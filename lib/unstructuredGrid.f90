module unstructuredGrid_m
    use unstructuredElement_m
    use kdTree_m
    use vector_m
    implicit none
    private

    type, extends(cell_t) :: cell_inFlow_t
        !! 流れ場格子構造体

        integer, allocatable :: boundFaceID(:)
            !! 境界面ID配列。内部格子なら境界面を持たず、配列サイズはゼロ。

        integer, allocatable :: adjacentCellID(:)
            !! 隣接セルID配列。

        real center(3)
            !! セル重心

        real flowVelocity(3)
            !! 流速

        real threshold
            !! 近傍セル判定の閾値

        integer, allocatable :: faceID(:)

    end type

    type boundaryTriangle_t
        !! 境界三角形面構造体

        integer nodeID(3)
            !! 節点ID配列

        real center(3)
            !! 面重心

        real normalVector(3)
            !! 法線ベクトル

        real moveVector(3)
            !! 移動量ベクトル（格子移動がある場合に必要になる）

    end type

    type, public :: FlowFieldUnstructuredGrid
        !! 流れ場非構造格子クラス
        private

        type(node_t), allocatable :: NODEs(:)
            !! 節点配列

        type(cell_inFlow_t), allocatable :: CELLs(:)
            !! セル配列

        type(boundaryTriangle_t), allocatable :: BoundFACEs(:)
            !! 境界面配列

        type(boundaryTriangle_t), allocatable :: FACEs(:)

        type(kdTree) kd_tree
            !! kd-tree（近傍セル探索用）

        real MIN_CDN(3)
            !! 座標の最小値(xyz)
        real MAX_CDN(3)
            !! 座標の最大値(xyz)

        integer :: num_refCellSearchFalse = 0
            !! 参照セル探索結果が悪いと判断された回数
        integer :: num_refCellSearch = 0
            !! 参照セル探索が行われた回数

        logical :: fph_flag = .false.

        contains
        private

        procedure, public :: nearest_cell, nearcell_check
        procedure, public :: nearest_search_exact, nearest_search_kdTree
        procedure, public :: get_flowVelocityInCELL, get_movementVectorOfBoundarySurface
        procedure, public :: get_cellCenterOf, get_allOfCellCenters
        procedure, public :: get_cellVerticesOf, get_MinMaxOfGrid
        procedure, public :: get_info => get_gridInformation

        procedure set_cellCenter, set_cellThreshold, set_MinMaxCDN, point2cellVelocity
        procedure, public :: read_VTK, read_array, read_INP, read_FLD, read_FPH

        !=====================================================================

        procedure, public :: setupWithFlowFieldFile, updateWithFlowFieldFile
        procedure nearer_cell
        procedure, public :: search_refCELL
        procedure, public :: adhesionCheckOnBound
        procedure, public ::  get_num_nearerSearchFalse, get_nearerSearchFalseRate

        procedure AdjacencySolvingProcess
        procedure read_adjacency, read_boundaries, solve_adjacencyOnFlowFieldUnstructuredGrid
        procedure output_boundaries, output_adjacency, boundary_setting, output_STL
        procedure read_cell2face
        
        ! procedure setup_kdTree

    end type

    public FlowFieldUnstructuredGrid_, FlowFieldUnstructuredGrid_withMeshFile

    contains

    function FlowFieldUnstructuredGrid_(FlowFieldFile) result(grid)
        !! 流れ場のコンストラクタ
        !! 流れ場ファイルの読み込み、前処理、kd-treeの構築を行う
        use path_operator_m
        character(*), intent(in) :: FlowFieldFile
            !! 流れ場ファイル名
        character(:), allocatable :: Dir
        type(FlowFieldUnstructuredGrid) grid

        call grid%setupWithFlowFieldFile(FlowFieldFile)
        call get_DirFromPath(FlowFieldFile, Dir)
        call grid%AdjacencySolvingProcess(Dir)    !流れ場の前処理
        ! call grid%setup_kdTree(Dir)

    end function

    function FlowFieldUnstructuredGrid_withMeshFile(FlowFieldFile, meshFile) result(grid)
        !! 流れ場のコンストラクタ(meshファイルと流速データファイルが分かれている場合)
        !! 流れ場ファイルの読み込み、前処理、kd-treeの構築を行う
        use path_operator_m
        character(*), intent(in) :: FlowFieldFile
            !! 流速データファイル名
        character(*), intent(in) :: meshFile
            !! メッシュファイル
        character(:), allocatable :: Dir
        type(FlowFieldUnstructuredGrid) grid

        call grid%setupWithFlowFieldFile(FlowFieldFile, meshFile)
        call get_DirFromPath(meshFile, Dir)
        call grid%AdjacencySolvingProcess(Dir)    !流れ場の前処理
        ! call grid%setup_kdTree(Dir)

    end function

    subroutine setupWithFlowFieldFile(self, FNAME, meshFile)
        !! 流れ場ファイルの読み込み
        !! VTK, INP, FLDに対応
        !! 独自フォーマットのArrayにも対応
        !! セル重心の算出も行う
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: FNAME
        character(*), intent(in), optional :: meshFile

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME, meshOnly=.false.)

            case('inp')
                call self%read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .true., findVelocity = .true.)

                if(.not.allocated(self%CELLs)) then !まだ未割り当てのとき
                    block
                        character(:), allocatable :: topologyFNAME
                        topologyFNAME = FNAME( : index(FNAME, '_', back=.true.)) // '0' // '.fld' !ゼロ番にアクセス
                        call self%read_FLD(topologyFNAME, findTopology= .true., findVelocity = .false.)
                    end block

                    call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)
                end if

            case('fph')
                ! call self%read_FPH(FNAME, findTopology= .true., findVelocity = .true.)

                if(.not.allocated(self%CELLs)) then !まだ未割り当てのとき
                    block
                        character(:), allocatable :: topologyFNAME
                        topologyFNAME = FNAME( : index(FNAME, '_', back=.true.)) // '0' // '.fph' !ゼロ番にアクセス
                        call self%read_FPH(topologyFNAME, findTopology= .true., findVelocity = .false.)
                    end block

                    call self%read_FPH(FNAME, findTopology= .false., findVelocity = .true.)
                end if

            case('array')
                call self%read_VTK(meshFile, meshOnly=.true.)
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', FNAME
                error stop
                    
        end select

        call self%set_MinMaxCDN()
        if(.not. self%fph_flag) call self%set_cellCenter()
        call self%set_cellThreshold()

    end subroutine

    subroutine updateWithFlowFieldFile(self, FNAME)
        !! 流れ場ファイルを読み込み、流れ場を更新する
        !! あくまで既存の流れ場の更新目的であり、セル数の異なるメッシュは想定していない
        !! セル重心の算出も行う
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: FNAME

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME, meshOnly=.false.)

            case('inp')
                call self%read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)

            case('fph')
                call self%read_FPH(FNAME, findTopology= .false., findVelocity = .true.)

            case('array')
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', FNAME
                error stop
                    
        end select

        call self%set_MinMaxCDN()
        if(.not. self%fph_flag) call self%set_cellCenter()
        call self%set_cellThreshold()

        call self%boundary_setting(first=.false.)
            
    end subroutine

    subroutine AdjacencySolvingProcess(self, dir)
        !! 前処理
        !! 隣接関係および境界面情報を解決する
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: dir
        logical success

        call self%read_adjacency(dir, success)
        if(success) then
            call self%read_boundaries(dir)

        else
            call self%solve_adjacencyOnFlowFieldUnstructuredGrid()
            call self%output_boundaries(dir)
            call self%output_adjacency(dir)

        end if

        call self%boundary_setting(first=.true.)

        if(.not. self%fph_flag) call self%output_STL(dir//'shape.stl')

    end subroutine

    subroutine set_MinMaxCDN(self)
        !! 節点群の座標最大最小の算出
        class(FlowFieldUnstructuredGrid) self
        real MinMax(3,2)

        MinMax = get_MinMaxCDN(self%NODEs)

        self%MIN_CDN = MinMax(:,1)
        self%MAX_CDN = MinMax(:,2)

    end subroutine

    subroutine get_MinMaxOfGrid(self, MIN_CDN, MAX_CDN)
        !! 節点群の座標最大最小を返す
        class(FlowFieldUnstructuredGrid) self
        real, intent(out) :: MIN_CDN(3), MAX_CDN(3)

        MIN_CDN = self%MIN_CDN
        MAX_CDN = self%MAX_CDN

    end subroutine

    subroutine set_cellCenter(self)
        !! セル重心の算出
        class(FlowFieldUnstructuredGrid) self
        integer II,IIMX, n, num_node, nodeID
        real vector(3)

        IIMX = size(self%CELLs)
        DO II = 1, IIMX
            num_node = size(self%CELLs(II)%nodeID)
            vector(:) = 0.0
            do n = 1, num_node
                nodeID = self%CELLs(II)%nodeID(n)
                vector(:) = vector(:) + self%NODEs(nodeID)%coordinate(:)
            end do
            self%CELLs(II)%center(:) = vector(:) / real(num_node)
        END DO

    end subroutine

    subroutine set_cellThreshold(self)
        !! セル閾値の算出
        class(FlowFieldUnstructuredGrid) self
        integer II,IIMX, n, num_node, nodeID, num_face, faceID
        real x, vector(3)

        IIMX = size(self%CELLs)

        if(.not. self%fph_flag) then

            DO II = 1, IIMX
                num_node = size(self%CELLs(II)%nodeID)
                x = 0.0
                do n = 1, num_node
                    nodeID = self%CELLs(II)%nodeID(n)
                    vector(:) = self%NODEs(nodeID)%coordinate(:) - self%CELLs(II)%center(:)
                    x = max(x, norm2(vector)) 
                end do
                self%CELLs(II)%threshold = x
                ! print*, self%CELLs(II)%threshold
            END DO

        else

            do II = 1, IIMX
                num_face = size(self%CELLs(II)%faceID)
                x = 0.0
                do n = 1, num_face
                    faceID = self%CELLs(II)%faceID(n)
                    vector(:) = self%FACEs(faceID)%center(:) - self%CELLs(II)%center(:)
                    x = max(x, norm2(vector)) 
                end do
                self%CELLs(II)%threshold = x

            end do

        end if

    end subroutine

    subroutine read_VTK(self, FNAME, meshOnly)
        !! VTKファイルから流れ場を取得する
        use VTK_operator_m
        class(FlowFieldUnstructuredGrid) self
        type(UnstructuredGrid_inVTK) vtk_mesh
        character(*), intent(in) :: FNAME
            !! ファイル名

        logical, intent(in) :: meshOnly
            !! メッシュだけ読み込み、流速などは無視するフラグ

        real, allocatable :: velocity(:,:)
        integer II, IIMX

        if(meshOnly) then
            call vtk_mesh%read(FNAME)
        else
            call vtk_mesh%read(FNAME, cellVector=velocity)
        end if

        self%NODEs = vtk_mesh%node_array
        
        IIMX = size(vtk_mesh%cell_array)
        if(.not.allocated(self%CELLs)) allocate(self%CELLs(IIMX))
        do II = 1, IIMX
            self%CELLs(II)%nodeID = vtk_mesh%cell_array(II)%nodeID
        end do

        if(.not.meshOnly) then
            do II = 1, IIMX
                self%CELLs(II)%flowVelocity(:) = velocity(:,II)
            end do
        endif

        ! print*, NODEs(KKMX)%coordinate(:)
        ! print*, CELLs(IIMX)%flowVelocity(:)
            
    end subroutine

    subroutine read_Array(self, FNAME)
        !! 独自フォーマットArrayファイルから、流速を読み込む
        use array_m
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: FNAME
        real, allocatable :: velocity(:,:)
        integer II

        call read_2dArray_asBinary(FNAME, velocity)

        if(size(self%CELLS) /= size(velocity, dim=2)) then
            print*, 'SIZE ERROR:', size(self%CELLS), size(velocity, dim=2)
            error stop
        end if

        do II = 1, size(self%CELLS)
            self%CELLs(II)%flowVelocity(:) = velocity(:,II)
            ! print*, velocity(:,II)
        end do

        ! print*, NODEs(KKMX)%coordinate(:)
        ! print*, CELLs(IIMX)%flowVelocity(:)
            
    end subroutine

    subroutine read_INP(self, FNAME)
        !! INPファイルを読み込み、節点データを要素データに変換する
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: FNAME
        INTEGER II,II2,KK,AAmax, IIMX2, AA, n_unit
        integer KKMX, IIMX
        character(6) cellshape
        real, allocatable :: UVWK(:,:)
        integer, allocatable :: ICN2(:,:)
        integer, allocatable :: CELL_TYPE2(:)

        print*, 'READ_INP:', trim(FNAME)
            
        open(newunit=n_unit,FILE=FNAME,status='old', action='read')
            read(n_unit,*)KKMX,IIMX2
            print*,'KKMX,IIMX2=',KKMX,IIMX2
            
            if(.not.allocated(self%NODEs)) allocate(self%NODEs(KKMX))
            allocate(ICN2(6,IIMX2), source=0)
            allocate(CELL_TYPE2(IIMX2), source=0)
                
            DO KK = 1, KKMX
                read(n_unit,*)AA, self%NODEs(KK)%coordinate(:)
            END DO
                
            II = 0
            DO II2 = 1, IIMX2
    
                read(n_unit,fmt='(I10)',advance='no')AA  !ここはセル番号なので無視
                read(n_unit,fmt='(I6)',advance='no')AA  !ここはなんかしらん（だいたいゼロ）
                read(n_unit,fmt='(A6)',advance='no')cellshape
    
                cellshape = adjustl(cellshape)    !左詰め
                IF ((cellshape=='tet').or.(cellshape=='prism').or.(cellshape=='pyr')) THEN
                    II = II +1
        
                    if (cellshape=='tet') then
                        CELL_TYPE2(II) = 0
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II)    
                        if (ICN2(1,II)==0.or.ICN2(4,II)==0) print*, 'ICN2_WARNING_tet:', ICN2(:,II)
                    
                    ELSE IF(cellshape=='prism') THEN
                        CELL_TYPE2(II) = 1
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II),ICN2(5,II),ICN2(6,II)
                        if (ICN2(1,II)==0.or.ICN2(6,II)==0) print*, 'ICN2_WARNING_prism:', ICN2(:,II)
        
                    ELSE IF(cellshape=='pyr') THEN
                        CELL_TYPE2(II) = 2
                        read(n_unit,*)ICN2(5,II),ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II) !INPは最初が山頂点であり、VTKでは最後が山頂点のため、読み込む順がこうなる。
                        if (ICN2(1,II)==0.or.ICN2(5,II)==0) print*, 'ICN2_WARNING_pyr:', ICN2(:,II)
        
                    end if
    
                ELSE
                    read(n_unit,'()')  !テトラでもプリズムでもピラミッドでもないならスルー
    
                ENDIF
    
            END DO

            IIMX = II
        
            allocate(UVWK(3,KKMX))
    
            read(n_unit,*)AAmax
            read(n_unit,'()')
            DO II = 1,AAmax
                read(n_unit,'()')
            END DO
            DO KK = 1, KKMX
                read(n_unit,*)AA, UVWK(:,KK)
            END DO
                
        close(n_unit)

        if(.not.allocated(self%CELLs)) then
            allocate(self%CELLs(IIMX))
            do II = 1, IIMX
                select case(CELL_TYPE2(II))
                    case(0)
                        self%CELLs(II)%nodeID = ICN2(1:4, II)
                    case(1)
                        self%CELLs(II)%nodeID = ICN2(1:6, II)
                    case(2)
                        self%CELLs(II)%nodeID = ICN2(1:5, II)
                end select
            end do
        end if
            
        call self%point2cellVelocity(UVWK)
            
    end subroutine

    subroutine read_FLD(self, FNAME, findTopology, findVelocity)
        !! FLDファイルから流れ場を取得する
        use SCT_file_reader_m
        class(FlowFieldUnstructuredGrid) self
        type(sct_grid_t) grid
        integer ii, iitet, iiwed, iipyr, iihex, iimx, iicnt
        integer kk, kkmx
        integer,allocatable :: tetras(:,:), wedges(:,:), pyramids(:,:), hexas(:,:)
        real(8),allocatable :: points(:,:)
        real(8),allocatable :: velocity(:,:)!, pressure(:)

        character(*), intent(in) :: FNAME
            !! ファイル名

        logical, intent(in) :: findTopology
            !! トポロジー情報を取得するフラグ

        logical, intent(in) :: findVelocity
            !! 流速情報を取得するフラグ

        print*, 'readFLD : ', trim(FNAME)

        call grid%read_SCT_file(FNAME)
        
        if(findTopology) then
            !ファイルが存在し, かつトポロジー情報が存在する場合以下の処理が行われる.  
            call grid%extract_cell_vertices(tetras, pyramids, wedges, hexas)
            if(allocated(tetras)) then
                call grid%get_2d_array_of_point_coords(points)
                iitet = grid%get_tetrahedron_count()
                iipyr = grid%get_pyramid_count()
                iiwed = grid%get_wedge_count()
                iihex = grid%get_hexahedron_count()
                iimx = grid%get_element_count()
                kkmx = grid%get_vertex_count()
    
                if(iihex>0) then
                    print*, 'Hexahedron is not yet supported.', iihex
                    error stop
                end if
    
                allocate(self%CELLs(iimx))
                allocate(self%NODEs(kkmx))
    
                do kk = 1, kkmx
                    self%NODEs(kk)%coordinate(:) = real(points(:,kk))
                end do
    
                iicnt = 1
                do ii = 1, iitet
                    self%CELLs(iicnt)%nodeID = tetras(:,ii)
                    iicnt = iicnt + 1
                end do
                do ii = 1, iiwed
                    self%CELLs(iicnt)%nodeID = wedges(:,ii)
                    iicnt = iicnt + 1
                end do
                do ii = 1, iipyr
                    self%CELLs(iicnt)%nodeID = pyramids(:,ii)
                    iicnt = iicnt + 1
                end do
    
            end if

        end if
        
        ! call grid%search_scalar_data("PRES",pressure)
        if(findVelocity) then
            call grid%search_vector_data("VEL",velocity)
            call self%point2cellVelocity(real(velocity))
        end if
        
    end subroutine

    subroutine read_FPH(self, FNAME, findTopology, findVelocity)
        !! FPHファイルから流れ場を取得する
        use SCF_file_reader_m
        use path_operator_m
        class(FlowFieldUnstructuredGrid) self
        type(scf_grid_t) grid
        character(:),allocatable :: dir
        real(4),allocatable :: points(:,:), velocity(:,:), bound_center(:,:), face_center(:,:)
        logical is_adjacencyFile, is_cell2faceFile
        integer iimx, kkmx, jjmx, ii, jj, kk, JB, num_boundFaces, n_unit

        character(*), intent(in) :: FNAME
            !! ファイル名

        logical, intent(in) :: findTopology
            !! トポロジー情報を取得するフラグ

        logical, intent(in) :: findVelocity
            !! 流速情報を取得するフラグ

        self%fph_flag = .true.

        print*, 'readFPH : ', trim(FNAME)

        call grid%read_SCF_file(FNAME)

        call get_DirFromPath(FNAME,dir)

        if(findTopology) then
            iimx = grid%get_fph_element_count()
            kkmx = grid%get_fph_vertex_count()
            jjmx = grid%get_fph_face_count()

            call grid%get_fph_2d_array_of_point_coords(points)
            call grid%get_face2vertices()
            call grid%get_face2cells()

            call grid%get_fph_faceCenter(face_center)
            call grid%get_fph_boundFaceIDs(num_boundFaces)
            call grid%get_fph_boundFaceCenter(bound_center)
            
            allocate(self%CELLs(iimx))
            allocate(self%NODEs(kkmx))
            allocate(self%FACEs(jjmx))
            allocate(self%BoundFACEs(num_boundFaces))

            inquire(file = dir//'adjacency.txt', exist = is_adjacencyFile)
            inquire(file = dir//'cell2face.txt', exist = is_cell2faceFile)

            do kk = 1, kkmx
                self%NODEs(kk)%coordinate(:) = real(points(:,kk))
            end do

            do jj = 1, jjmx
                self%FACEs(jj)%center(:) = real(face_center(:,jj))
            end do

            do JB = 1, num_boundFaces
                self%BoundFACEs(JB)%center(:) = real(bound_center(:,JB))
            end do

            if(.not. is_cell2faceFile) then

                call grid%get_cell2faces()                
                call grid%output_fph_cell2face(dir)

            end if

            call self%read_cell2face(dir)

            if(.not. is_adjacencyFile) then

                call grid%get_cell_offsets()
                call grid%get_cell2boundFace()
                call grid%get_fph_adjacentCellIDs()

                call grid%output_fph_boundFace(dir)
                call grid%output_fph_adjacentCell(dir)
                call grid%output_fph_vtk(dir)

            end if

        end if

        if(findVelocity) then

            call grid%search_fph_vector_data("VEL",velocity)

            iimx = size(self%CELLs)            
            do ii = 1, iimx
                self%CELLs(II)%flowVelocity(:) = real(velocity(:,II))
            end do

        end if
        
    end subroutine

    subroutine read_cell2face(self,dir)
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: dir
        integer n_unit, num_cell2face, NF, ii, iimx
        character(255) str

        iimx = size(self%CELLs)

        open(newunit = n_unit, file = dir//'cell2face.txt', status = 'old')
            do ii = 1, iimx
                read(n_unit,'(A)') str
                read(str, *) num_cell2face
                allocate(self%CELLs(II)%faceID(num_cell2face))
                read(str, *) NF, self%CELLs(II)%faceID(:)
            end do
        close(n_unit)

    end subroutine

    subroutine read_adjacency(self, path, success)
        !! セルの隣接関係情報をTXTファイルから読み込む
        use filename_m, only : adjacencyFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
            !! 隣接関係ファイルへのパス

        logical, intent(out) :: success
            !! 指定ファイルから正しく情報を取得できたら`.true.`が返る

        integer II,NA, n_unit, num_cells, num_adj, num_BF, NCMAX
        character(:), allocatable :: FNAME
        character(255) str
                
        FNAME = trim(path)//adjacencyFileName
        inquire(file = FNAME, exist = success)
        if(.not.success) then
            print*, 'AdjacencyFile was not found.'
            return
        end if

        print*, 'READ : ', FNAME

        open(newunit=n_unit, FILE=FNAME, status='old', action='read')
            read(n_unit,*) num_cells

            if(num_cells /= size(self%CELLs)) then
                print*, '**SIZE MISMATCH** :', num_cells, size(self%CELLs)
                success = .false.
                return
            end if

            read(n_unit,*) NCMAX

            DO II = 1, num_cells
                read(n_unit,'(A)') str
                read(str, *) num_adj
                allocate(self%CELLs(II)%adjacentCellID(num_adj))
                read(str, *) NA, self%CELLs(II)%adjacentCellID(:)
            END DO

            DO II = 1, num_cells
                read(n_unit,'(A)') str
                read(str, *) num_BF
                allocate(self%CELLs(II)%boundFaceID(num_BF))
                read(str, *) NA, self%CELLs(II)%boundFaceID(:)
            END DO

        close(n_unit)

    end subroutine

    subroutine output_adjacency(self, path)
        !! セルの隣接関係情報をTXTファイルに出力
        use filename_m, only : adjacencyFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        integer II, n_unit, num_cells, NCMAX
        character(:), allocatable :: FNAME
                
        FNAME = trim(path)//adjacencyFileName

        print*, 'OUTPUT:', FNAME

        num_cells = size(self%CELLs)
        NCMAX = 5   !size(NEXT_CELL(:,:), dim=1)

        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit,*) num_cells
            write(n_unit,*) NCMAX

            DO II = 1, num_cells
                write(n_unit,'(*(i0:,X))') size(self%CELLs(II)%adjacentCellID), self%CELLs(II)%adjacentCellID(:)
            END DO

            DO II = 1, num_cells
                write(n_unit,'(*(i0:,X))') size(self%CELLs(II)%boundFaceID), self%CELLs(II)%boundFaceID(:)
            END DO

        close(n_unit)

    end subroutine

    subroutine read_boundaries(self, path)
        !! 境界面情報をTXTファイルから読み込む
        use filename_m, only : boundaryFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(:), allocatable :: FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'READ : ', FNAME
        open(newunit=n_unit, FILE=FNAME , status='old', action='read')
            read(n_unit,*) JBMX
            if(.not. allocated(self%BoundFACEs)) allocate(self%BoundFACEs(JBMX))
            do JB = 1, JBMX
                read(n_unit,*) self%BoundFACEs(JB)%nodeID(:)
            end do
        close(n_unit)
        
    end subroutine

    subroutine output_boundaries(self, path)
        !! 境界面情報をTXTファイルに出力
        use filename_m, only : boundaryFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(:), allocatable :: FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'OUTPUT:', FNAME
        JBMX = size(self%BoundFACEs)
        open(newunit=n_unit, FILE=FNAME , STATUS='replace')
            write(n_unit,*) JBMX
            do JB = 1, JBMX
                write(n_unit,'(*(i0:,X))') self%BoundFACEs(JB)%nodeID(:)
            end do
        close(n_unit)
        
    end subroutine

    subroutine setup_kdTree(self, path)
        !! kd-treeの構築
        use filename_m, only : kdTreeFName 
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        character(:), allocatable :: FNAME
        real, allocatable :: xyz(:,:)
        logical existance

        FNAME = trim(path)//kdTreeFName

        xyz = self%get_allOfCellCenters()

        inquire(file = FNAME, exist=existance)
        if(.not.existance) then

            self%kd_tree = kdTree_(xyz)
            ! 時間かかるので、今は使わない
            ! call self%kd_tree%saveAsTXT(FNAME)
            ! print*, 'OUTPUT kdtree:', FNAME

        else

            call self%kd_tree%read_kdTree(FNAME)
            print*, 'READ kdtree:', FNAME

        end if

    end subroutine

    function nearest_cell(self, X) result(nearestCellID)
        !! 最近傍セル探索
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        real, intent(in) :: X(3)
            !! 探索対象座標
        integer nearestCellID

        !!@note 厳密探索かkdツリー探索かはここで切り替える

        nearestCellID = self%nearest_search_exact(X)
        ! nearestCellID = self%nearest_search_kdTree(X)

    end function

    function nearest_search_exact(self, X) result(nearestCellID)
        !! 厳密最近傍セル探索
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        real, intent(in) :: X(3)
            !! 探索対象座標
        integer II, IIMX
        real, allocatable :: distance(:)
        integer nearestCellID

        IIMX = size(self%CELLs)
        allocate(distance(IIMX))

        !$omp parallel do
        DO II = 1,IIMX
            ! distance(II) = norm2(self%CELLs(II)%center(:) - X(:))
            distance(II) = norm2_squared(self%CELLs(II)%center(:) - X(:))
        END DO
        !$omp end parallel do 
        
        nearestCellID = minloc(distance, dim=1)   !最小値インデックス

    end function

    function nearest_search_kdTree(self, X) result(nearestCellID)
        !!kdツリーによる最近傍セル探索
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        real, intent(in) :: X(3)
            !! 探索対象座標
        real, allocatable :: xyz(:,:)
        integer nearestCellID

        xyz = self%get_allOfCellCenters()

        call self%kd_tree%search(xyz, X, nearestCellID)
        
    end function

    function nearer_cell(self, X, NCN) result(nearCellID)
        !! 隣接セルを起点に近傍セル探索
        !! セルの隣接関係が悪いと上手く探索できないおそれあり
        class(FlowFieldUnstructuredGrid) self
        integer, intent(in) :: NCN
            !! 探索開始セルID
        real, intent(in) :: X(3)
            !! 探索対象座標
        integer NA, featuringCellID, adjaCellID
        integer, allocatable :: adjacentCellIDs(:)
        real distance, distance_min
        logical update
        integer nearCellID

        nearCellID = NCN
        distance_min = norm2(self%CELLs(nearCellID)%center(:) - X(:))   !注目セル重心と粒子との距離
        update = .true.

        do while(update)    !更新が起こり続ける限り繰り返し
            update = .false.

            featuringCellID = nearCellID       !注目セル

            adjacentCellIDs = self%CELLs(featuringCellID)%adjacentCellID(:)  !注目セルの全隣接セル

            do NA = 1, size(adjacentCellIDs)  !注目セルの全隣接セルに対してループ。

                adjaCellID = adjacentCellIDs(NA)  !注目セルに隣接するセルのひとつ
                ! if (adjacentCELL <= 0) cycle checkAdjacent

                distance = norm2(self%CELLs(adjaCellID)%center(:) - X(:))   !隣接セル重心と粒子との距離
                if(distance < distance_min) then
                    nearCellID = adjaCellID
                    distance_min = distance
                    update = .true.
                end if

            end do

        end do
        
    end function

    function nearcell_check(self, X, NCN) result(isNear)
        !! 近傍セル探索の結果が妥当かどうかをチェック
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: X(3)
            !! 探索対象座標
        integer, intent(in) :: NCN
            !! 近傍探索結果セルID
        real :: distance
        logical isNear

        distance = norm2(X(:) - self%CELLs(NCN)%center(:))

        !遠くのセルを参照していないかどうかのチェック
        !参照セルとの距離がセル閾値未満であればOK（この条件は経験則でしかない）
        isNear = (distance < self%CELLs(NCN)%threshold)

    end function

    subroutine search_refCELL(self, X, reference_cell, stat)
        !! 参照セル探索
        !! 主に近傍探索が呼ばれるが、探索が芳しくない場合は最近傍探索が呼ばれる
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: X(3)
            !! 探索対象座標
        integer, intent(inout) :: reference_cell
            !! 参照セル
        logical, optional :: stat
            !! 

        self%num_refCellSearch = self%num_refCellSearch + 1

        reference_cell = self%nearer_cell(X, reference_cell)
        if(present(stat)) stat = .True.

        if (.not.self%nearcell_check(X(:), reference_cell)) then
            reference_cell = self%nearest_cell(X)
            if(present(stat)) stat = .false.
            self%num_refCellSearchFalse = self%num_refCellSearchFalse + 1
        end if
    
    end subroutine

    function get_num_nearerSearchFalse(self) result(num_nearerSearchFalse)
        !! 近傍セル探索の結果が悪いと判断された回数を返す
        class(FlowFieldUnstructuredGrid) self
        integer num_nearerSearchFalse

        num_nearerSearchFalse = self%num_refCellSearchFalse

    end function

    function get_nearerSearchFalseRate(self) result(num_nearerSearchFalseRate)
        !! 近傍セル探索の結果が悪いと判断された比率を返す
        class(FlowFieldUnstructuredGrid) self
        real num_nearerSearchFalseRate

        num_nearerSearchFalseRate = 100. * real(self%num_refCellSearchFalse) / real(self%num_refCellSearch + 1)

    end function
                     
    subroutine boundary_setting(self, first)
        !! 全境界面に対して外向き法線ベクトルと重心を算出
        class(FlowFieldUnstructuredGrid) self
        logical, intent(in) :: first
            !! 初期ステップであるか否か
        integer II, JJ, JB, IIMX, JBMX, nodeID(3)
        real :: a(3), b(3), r(3), normalVector(3)
        type(boundaryTriangle_t), allocatable :: BoundFACEs_pre(:)

        IIMX = size(self%CELLs)

        if(.not.first) BoundFACEs_pre = self%BoundFACEs
        
        do II = 1, IIMX
            
            do JJ = 1, size(self%CELLs(II)%boundFaceID)
                JB = self%CELLs(II)%boundFaceID(JJ)
                nodeID(:) = self%BoundFACEs(JB)%nodeID(:)
                
                if(.not. self%fph_flag) then
                    self%BoundFACEs(JB)%center(:) = ( self%NODEs(nodeID(1))%coordinate(:) &
                                                + self%NODEs(nodeID(2))%coordinate(:) &
                                                + self%NODEs(nodeID(3))%coordinate(:) ) / 3.0
                end if

                a(:) =  self%NODEs(nodeID(2))%coordinate(:) - self%NODEs(nodeID(1))%coordinate(:)
                b(:) =  self%NODEs(nodeID(3))%coordinate(:) - self%NODEs(nodeID(1))%coordinate(:)
                normalVector(:) = cross_product(a, b)

                normalVector(:) = normalize_vector(normalVector(:))
            
                r(:) = self%CELLs(II)%center(:) - self%BoundFACEs(JB)%center(:)  !面重心からセル重心へのベクトル
                if(dot_product(normalVector(:), r(:)) > 0.0) then
                    normalVector(:) = normalVector(:) * (-1.0) !内積が正なら内向きなので、外に向ける
                end if

                self%BoundFACEs(JB)%normalVector(:) = normalVector(:)
                ! print*,'center:',BoundFACEs(JB)%center(:)
                ! print*,'n_vector:',BoundFACEs(JB)%normalVector(:)
            end do
        
        end do

        ! movevectorが上手くいってない?
        JBMX = size(self%BoundFACEs)
        if(first) then
            do JB = 1, JBMX
                self%BoundFACEs(JB)%moveVector(:) = 0.0
            end do
        else
            do JB = 1, JBMX
                self%BoundFACEs(JB)%moveVector(:) = self%BoundFACEs(JB)%center(:) - BoundFACEs_pre(JB)%center(:)
            end do
        end if

    end subroutine

    subroutine adhesionCheckOnBound(self, position, radius, cellID, stat)
        !! 境界面への飛沫付着判定
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        double precision, intent(in) :: position(3)
            !! 飛沫座標
        double precision, intent(in) :: radius
            !! 飛沫半径
        integer, intent(in) :: cellID
            !! 判定対象セルID
        integer, intent(out) :: stat
            !! 付着が起こらなければゼロ、起これば付着面の境界面IDが返る
        integer JJ, JB

        double precision :: r_vector(3), inner

        stat = 0

        do JJ = 1, size(self%CELLs(CellID)%boundFaceID)
            JB = self%CELLs(CellID)%boundFaceID(JJ)

            r_vector(:) = position(:) - self%BoundFACEs(JB)%center(:)

            inner = dot_product(r_vector(:), self%BoundFACEs(JB)%normalVector(:))
            !外向き法線ベクトルと位置ベクトルの内積は、平面からの飛び出し量に相当

            if (inner + radius > 0.d0) then !(飛び出し量+飛沫半径)がゼロ以上なら付着判定  
                stat = JB       !付着した境界面番号
            end if
        end do

    end subroutine

    subroutine point2cellVelocity(self, pointVector)
        !! 節点定義のベクトル場を、セル重心定義に換算
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: pointVector(:,:)
            !! 節点定義のベクトル配列（ベクトル長さ x 節点数）
        integer II, IIMX, n, ID, num_node

        if(.not.allocated(self%CELLs)) then
            print*, '**MISS point2cellVelocity** CELL_ARRAY is not yet allocated.'
            return
        end if

        IIMX = size(self%CELLs)
        DO II = 1, IIMX
            self%CELLs(II)%flowVelocity(:) = 0.0
            num_node = size(self%CELLs(II)%nodeID)
            do n = 1, num_node
                ID = self%CELLs(II)%nodeID(n)
                self%CELLs(II)%flowVelocity(:) = self%CELLs(II)%flowVelocity(:) + pointVector(:,ID)
            end do
            self%CELLs(II)%flowVelocity(:) = self%CELLs(II)%flowVelocity(:) / real(num_node)
        END DO
  
    end subroutine

    subroutine output_STL(self, fname)
        !! 境界面形状情報をSTL形式で出力
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: fname
        integer n_unit, JB, nodeID(3)
        real cross(3)

        print*, 'output_STL : ', fname

        open(newunit=n_unit, file=fname, status='replace')
            write(n_unit, '("solid test")')

            do JB = 1, size(self%BoundFACEs)
                write(n_unit, '(" facet normal", 3(X,F7.4))') self%BoundFACEs(JB)%normalVector(:)

                nodeID = self%BoundFACEs(JB)%nodeID(:)

                block
                    real, dimension(3) :: a,b

                    a = self%NODEs(nodeID(2))%coordinate(:) - self%NODEs(nodeID(1))%coordinate(:)
                    b = self%NODEs(nodeID(3))%coordinate(:) - self%NODEs(nodeID(1))%coordinate(:)
                    cross = cross_product(a, b)
                end block

                write(n_unit, '(" outer loop")')

                !面の向きによって点の並びを変えないと表示が上手く行かない（？）
                if (dot_product(cross, self%BoundFACEs(JB)%normalVector(:)) > 0.) then
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(1))%coordinate(:)
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(2))%coordinate(:)
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(3))%coordinate(:)
                else
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(1))%coordinate(:)
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(3))%coordinate(:)
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(nodeID(2))%coordinate(:)
                end if

                write(n_unit, '(" endloop")')
                write(n_unit, '(" endfacet")')
            end do

            write(n_unit, '("endsolid test")')
        close(n_unit)

    end subroutine

    subroutine solve_adjacencyOnFlowFieldUnstructuredGrid(self)
        !!境界面と隣接関係を解決し、結果を非構造格子クラスに格納
        use adjacencySolver_m
        class(FlowFieldUnstructuredGrid) self
        integer i, j, num_adjacent, num_boundFace, num_node
        integer max_vertex
        integer, allocatable :: cellVertices(:,:)
        integer, allocatable :: adjacentCellIDArray(:,:)
        integer, allocatable :: cellBoundFaces(:,:)
        integer, allocatable :: triangleBoundFaceVertices(:,:)

        max_vertex = 0
        do i = 1, size(self%CELLs)
            num_node = size(self%CELLs(i)%nodeID(:))
            max_vertex = max(max_vertex, num_node)  !頂点数の最大値の探索
        end do
        allocate(cellVertices(max_vertex, size(self%CELLs)), source=None)
        do i = 1, size(self%CELLs)
            num_node = size(self%CELLs(i)%nodeID(:))
            cellVertices(1:num_node, i) = self%CELLs(i)%nodeID(1:num_node)
        end do

        call solve_BoundaryAndAdjacency(cellVertices, cellBoundFaces, triangleBoundFaceVertices, adjacentCellIDArray)

        allocate(self%BoundFACEs(size(triangleBoundFaceVertices, dim=2)))
        do j = 1, size(self%BoundFACEs)
            self%BoundFACEs(j)%nodeID = triangleBoundFaceVertices(:,j)
        end do

        do i = 1, size(self%CELLs)
            num_boundFace = count(cellBoundFaces(:,i)/=None)
            self%CELLs(i)%boundFaceID = cellBoundFaces(1:num_boundFace, i)

            num_adjacent = count(adjacentCellIDArray(:,i)/=None)
            self%CELLs(i)%adjacentCellID = adjacentCellIDArray(1:num_adjacent, i)
        end do

    end subroutine

    function extensionOf(FileName) result(extension)
        !! ファイル名の拡張子を返す
        character(*), intent(in) :: FileName
        character(:), allocatable :: extension

        extension = FileName(index(FileName, '.', back=.true.)+1 : )

    end function
    
    function get_flowVelocityInCELL(self, ID) result(velocity)
        !! 指定IDにおけるセルの流速を返す
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real velocity(3)

        velocity = self%CELLs(ID)%flowVelocity

    end function

    function get_movementVectorOfBoundarySurface(self, ID) result(vector)
        !! 指定IDにおける境界面の移動量ベクトルを返す
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real vector(3)

        vector = self%BoundFACEs(ID)%moveVector(:)

    end function

    function get_cellCenterOf(self, ID) result(center)
        !! 指定IDにおけるセル重心
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real center(3)

        center = self%CELLs(ID)%center

    end function

    function get_allOfCellCenters(self) result(centers)
        !! 全セルの重心座標を２次元配列で返す
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        real, allocatable :: centers(:,:)
        integer i, num_cell

        num_cell = size(self%CELLs)
        allocate(centers(3,num_cell))
        do i = 1, num_cell
            centers(:,i) = self%CELLs(i)%center
        end do
    end function

    function get_cellVerticesOf(self, ID) result(vertices)
        !! 指定IDのセルにおける頂点座標配列を返す
        !! ２次元配列：（xyz、頂点）
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real, allocatable :: vertices(:,:)
        integer i, num_node, nodeID

        num_node = size(self%CELLs(ID)%nodeID)
        allocate(vertices(3, num_node))
        do i = 1, num_node
            nodeID = self%CELLs(ID)%nodeID(i)
            vertices(:,i) = self%NODEs(nodeID)%coordinate
        end do

    end function

    function get_gridInformation(self, name) result(info)
        !! 節点数もしくはセル数を取得
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        character(*), intent(in) :: name
            !! 'node' or 'cell'
        integer info

        select case(name)
        case('node')
            info = size(self%NODEs)
        case('cell')
            info = size(self%CELLs)
        end select

    end function

    
end module unstructuredGrid_m