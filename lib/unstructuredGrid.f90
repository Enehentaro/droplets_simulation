module unstructuredGrid_m
    use unstructuredElement_m
    implicit none
    private

    type, extends(cell_t) :: cell_inFlow_t
        integer, allocatable :: boundFaceID(:), adjacentCellID(:)
        real center(3), flowVelocity(3), threshold
    end type

    type boundaryTriangle_t
        integer nodeID(3)
        real center(3), normalVector(3), moveVector(3)
    end type

    type, public :: FlowFieldUnstructuredGrid
        private
        type(node_t), allocatable :: NODEs(:)
        type(cell_inFlow_t), allocatable :: CELLs(:)
        type(boundaryTriangle_t), allocatable :: BoundFACEs(:)

        real MIN_CDN(3), MAX_CDN(3)
        integer :: num_refCellSearchFalse = 0, num_refCellSearch = 0

        contains
        private

        procedure, public :: nearest_cell, nearcell_check, get_MinMaxOfGrid
        procedure, public :: get_flowVelocityInCELL, get_movementVectorOfBoundarySurface
        procedure, public :: get_cellCenterOf, get_allOfCellCenters
        procedure, public :: get_cellVerticesOf
        procedure, public :: get_info => get_gridInformation

        procedure set_cellCenter, set_cellThreshold, set_MinMaxCDN, point2cellVelocity
        procedure, public :: read_VTK, read_array, read_INP, read_FLD

        !=====================================================================

        procedure, public :: setupWithFlowFieldFile, updateWithFlowFieldFile
        procedure nearer_cell
        procedure, public :: adhesionCheckOnBound
        procedure, public :: refCellSearchInfo, search_refCELL

        procedure AdjacencySolvingProcess
        procedure read_adjacency, read_boundaries, solve_adacencyOnFlowFieldUnstructuredGrid
        procedure output_boundaries, output_adjacency, boundary_setting, output_STL

    end type

    public FlowFieldUnstructuredGrid_

    contains

    type(FlowFieldUnstructuredGrid) function FlowFieldUnstructuredGrid_(FlowFieldFile, meshFile)
        use path_operator_m
        character(*), intent(in) :: FlowFieldFile
        character(*), intent(in), optional :: meshFile
        character(:), allocatable :: Dir

        if(present(meshFile)) then
            call FlowFieldUnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile, meshFile)
            call get_DirFromPath(meshFile, Dir)
        else
            call FlowFieldUnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile)
            call get_DirFromPath(FlowFieldFile, Dir)
        end if

        call FlowFieldUnstructuredGrid_%AdjacencySolvingProcess(Dir)    !流れ場の前処理

    end function

    subroutine setupWithFlowFieldFile(self, FNAME, meshFile)
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

            case('array')
                call self%read_VTK(meshFile, meshOnly=.true.)
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', FNAME
                error stop
                    
        end select

        call self%set_MinMaxCDN()
        call self%set_cellCenter()
        call self%set_cellThreshold()

    end subroutine

    subroutine updateWithFlowFieldFile(self, FNAME)
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: FNAME

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME, meshOnly=.false.)

            case('inp')
                call self%read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)

            case('array')
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', FNAME
                error stop
                    
        end select

        call self%set_MinMaxCDN()
        call self%set_cellCenter()
        call self%set_cellThreshold()

        call self%boundary_setting(first=.false.)
            
    end subroutine

    subroutine AdjacencySolvingProcess(self, dir)
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: dir
        logical success

        call self%read_adjacency(dir, success)
        if(success) then
            call self%read_boundaries(dir)

        else
            call self%solve_adacencyOnFlowFieldUnstructuredGrid()
            call self%output_boundaries(dir)
            call self%output_adjacency(dir)

        end if

        call self%boundary_setting(first=.true.)

        call self%output_STL(dir//'shape.stl')

    end subroutine

    subroutine set_MinMaxCDN(self) !節点群の座標最大最小の算出
        class(FlowFieldUnstructuredGrid) self
        real MinMax(3,2)

        MinMax = get_MinMaxCDN(self%NODEs)

        self%MIN_CDN = MinMax(:,1)
        self%MAX_CDN = MinMax(:,2)

    end subroutine

    subroutine get_MinMaxOfGrid(self, MIN_CDN, MAX_CDN) !節点群の座標最大最小を返す
        class(FlowFieldUnstructuredGrid) self
        real, intent(out) :: MIN_CDN(3), MAX_CDN(3)

        MIN_CDN = self%MIN_CDN
        MAX_CDN = self%MAX_CDN

    end subroutine

    subroutine set_cellCenter(self) !セル重心の算出
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

    subroutine set_cellThreshold(self) !セル閾値の算出
        class(FlowFieldUnstructuredGrid) self
        integer II,IIMX, n, num_node, nodeID
        real x, vector(3)

        IIMX = size(self%CELLs)
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

    end subroutine

    subroutine read_VTK(self, FNAME, meshOnly)
        use VTK_operator_m
        class(FlowFieldUnstructuredGrid) self
        type(UnstructuredGrid_inVTK) vtk_mesh
        character(*), intent(in) :: FNAME
        logical, intent(in) :: meshOnly
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
        !  INPファイルを読み込み、節点データを要素データに変換する
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
        use SCT_file_reader_m
        class(FlowFieldUnstructuredGrid) self
        type(sct_grid_t) grid
        integer ii, iitet, iiwed, iipyr, iihex, iimx, iicnt
        integer kk, kkmx
        integer,allocatable :: tetras(:,:), wedges(:,:), pyramids(:,:), hexas(:,:)
        logical, intent(in) :: findTopology, findVelocity
        real(8),allocatable :: points(:,:)
        real(8),allocatable :: velocity(:,:)!, pressure(:)
        character(*), intent(in) :: FNAME

        print*, 'readFLD : ', trim(FNAME)

        call grid%read_SCT_file(FNAME)
        
        if(findTopology) then
            !!ファイルが存在し, かつトポロジー情報が存在する場合以下の処理が行われる.  
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

    subroutine read_adjacency(self, path, success)
        use filename_m, only : adjacencyFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        logical, intent(out) :: success
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
        use filename_m, only : boundaryFileName
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(:), allocatable :: FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'READ : ', FNAME
        open(newunit=n_unit, FILE=FNAME , status='old', action='read')
            read(n_unit,*) JBMX
            allocate(self%BoundFACEs(JBMX))
            do JB = 1, JBMX
                read(n_unit,*) self%BoundFACEs(JB)%nodeID(:)
            end do
        close(n_unit)
        
    end subroutine

    subroutine output_boundaries(self, path)
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

    integer function nearest_cell(self, X) !最も近いセルNCNの探索
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: X(3)
        integer II, IIMX
        real, allocatable :: distance(:)

        IIMX = size(self%CELLs)
        allocate(distance(IIMX))

        !$omp parallel do
        DO II = 1,IIMX
            distance(II) = norm2(self%CELLs(II)%center(:) - X(:))
        END DO
        !$omp end parallel do 
        
        nearest_cell = minloc(distance, dim=1)   !最小値インデックス
        ! print*, distance(nearest_cell)
        
    end function

    integer function nearer_cell(self, X, NCN)  !近セルの探索（隣接セルから）
        class(FlowFieldUnstructuredGrid) self
        integer, intent(in) :: NCN
        real, intent(in) :: X(3)
        integer NA, featuringCellID, adjaCellID
        integer, allocatable :: adjacentCellIDs(:)
        real distance, distance_min
        logical update

        nearer_cell = NCN
        distance_min = norm2(self%CELLs(nearer_cell)%center(:) - X(:))   !注目セル重心と粒子との距離
        update = .true.

        do while(update)    !更新が起こり続ける限り繰り返し
            update = .false.

            featuringCellID = nearer_cell       !注目セル

            adjacentCellIDs = self%CELLs(featuringCellID)%adjacentCellID(:)  !注目セルの全隣接セル

            do NA = 1, size(adjacentCellIDs)  !注目セルの全隣接セルに対してループ。

                adjaCellID = adjacentCellIDs(NA)  !注目セルに隣接するセルのひとつ
                ! if (adjacentCELL <= 0) cycle checkAdjacent

                distance = norm2(self%CELLs(adjaCellID)%center(:) - X(:))   !隣接セル重心と粒子との距離
                if(distance < distance_min) then
                    nearer_cell = adjaCellID
                    distance_min = distance
                    update = .true.
                end if

            end do

        end do
        
    end function

    logical function nearcell_check(self, X, NCN)
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: X(3)
        integer, intent(in) :: NCN
        real :: distance

        distance = norm2(X(:) - self%CELLs(NCN)%center(:))

        !遠くのセルを参照していないかどうかのチェック
        !参照セルとの距離がセル閾値未満であればOK（この条件は経験則でしかない）
        if (distance < self%CELLs(NCN)%threshold) then
            nearcell_check = .True.
        else
            nearcell_check = .False.
            ! print*, 'nearcell_check:False', distance, self%CELLs(NCN)%threshold
        end if

    end function

    subroutine search_refCELL(self, X, reference_cell, stat)
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: X(3)
        integer, intent(inout) :: reference_cell
        logical, optional :: stat

        self%num_refCellSearch = self%num_refCellSearch + 1

        reference_cell = self%nearer_cell(X, reference_cell)
        if(present(stat)) stat = .True.

        if (.not.self%nearcell_check(X(:), reference_cell)) then
            reference_cell = self%nearest_cell(X)
            if(present(stat)) stat = .false.
            self%num_refCellSearchFalse = self%num_refCellSearchFalse + 1
        end if
    
    end subroutine

    integer function refCellSearchInfo(self, name)
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: name

        select case(name)
            case('NumFalse')
                refCellSearchInfo = self%num_refCellSearchFalse
            case('FalseRate')
                refCellSearchInfo = 100 * self%num_refCellSearchFalse / (self%num_refCellSearch + 1)
            case default
                print*, 'ERROR refCellSearchInfo : ', name
                error stop
        end select

    end function
                     
    subroutine boundary_setting(self, first) !全境界面に対して外向き法線ベクトルと重心を算出
        use vector_m
        class(FlowFieldUnstructuredGrid) self
        logical, intent(in) :: first
        integer II, JJ, JB, IIMX, JBMX, nodeID(3)
        real :: a(3), b(3), r(3), normalVector(3)
        type(boundaryTriangle_t), allocatable :: BoundFACEs_pre(:)

        IIMX = size(self%CELLs)

        if(.not.first) BoundFACEs_pre = self%BoundFACEs
        
        do II = 1, IIMX
            
            do JJ = 1, size(self%CELLs(II)%boundFaceID)
                JB = self%CELLs(II)%boundFaceID(JJ)
                nodeID(:) = self%BoundFACEs(JB)%nodeID(:)
                
                self%BoundFACEs(JB)%center(:) = ( self%NODEs(nodeID(1))%coordinate(:) &
                                            + self%NODEs(nodeID(2))%coordinate(:) &
                                            + self%NODEs(nodeID(3))%coordinate(:) ) / 3.0
            
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
        use vector_m
        class(FlowFieldUnstructuredGrid) self
        double precision, intent(in) :: position(3), radius
        integer, intent(in) :: cellID
        integer, intent(out) :: stat
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
        class(FlowFieldUnstructuredGrid) self
        real, intent(in) :: pointVector(:,:)
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
        class(FlowFieldUnstructuredGrid) self
        character(*), intent(in) :: fname
        integer i, n_unit, JB

        print*, 'output_STL : ', fname

        open(newunit=n_unit, file=fname, status='replace')
            write(n_unit, '("solid test")')

            do JB = 1, size(self%BoundFACEs)
                write(n_unit, '(" facet normal", 3(X,F7.4))') self%BoundFACEs(JB)%normalVector(:)
                write(n_unit, '(" outer loop")')

                do i = 1, 3
                    write(n_unit, '("  vertex", 3(X,E11.4))') self%NODEs(self%BoundFACEs(JB)%nodeID(i))%coordinate(:)
                end do

                write(n_unit, '(" endloop")')
                write(n_unit, '(" endfacet")')
            end do

            write(n_unit, '("endsolid test")')
        close(n_unit)

    end subroutine

    subroutine solve_adacencyOnFlowFieldUnstructuredGrid(self)
        use adjacencySolver_m
        class(FlowFieldUnstructuredGrid) self
        integer i, j, num_adjacent, num_boundFace
        integer, parameter :: max_vertex=6, max_adjacent=4, max_boundFace=4
        integer, allocatable :: cellVertices(:,:)
        integer, allocatable :: adjacentCellArray(:,:)
        integer, allocatable :: cellBoundFaces(:,:)
        integer, allocatable :: boundFaceVertices(:,:)

        allocate(cellVertices(max_vertex, size(self%CELLs)), source=0)
        do i = 1, size(self%CELLs)
            cellVertices(1:size(self%CELLs(i)%nodeID(:)), i) = self%CELLs(i)%nodeID(:)
        end do

        allocate(cellBoundFaces(max_adjacent, size(self%CELLs)), source=0)
        allocate(adjacentCellArray(max_boundFace, size(self%CELLs)), source=0)
        call solve_BoundaryAndAdjacency(cellVertices, cellBoundFaces, boundFaceVertices, adjacentCellArray)

        allocate(self%BoundFACEs(size(boundFaceVertices, dim=2)))
        do j = 1, size(self%BoundFACEs)
            self%BoundFACEs(j)%nodeID = boundFaceVertices(:,j)
        end do

        do i = 1, size(self%CELLs)
            num_boundFace = max_boundFace - count(cellBoundFaces(:,i)==0)
            self%CELLs(i)%boundFaceID = cellBoundFaces(1:num_boundFace, i)

            num_adjacent = max_adjacent - count(adjacentCellArray(:,i)==0)
            self%CELLs(i)%adjacentCellID = adjacentCellArray(1:num_adjacent, i)
        end do

    end subroutine

    function extensionOf(FileName) result(extension)
        character(*), intent(in) :: FileName
        character(:), allocatable :: extension

        extension = FileName(index(FileName, '.', back=.true.)+1 : )

    end function
    
    function get_flowVelocityInCELL(self, ID) result(velocity)
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real velocity(3)

        velocity = self%CELLs(ID)%flowVelocity

    end function

    function get_movementVectorOfBoundarySurface(self, ID) result(vector)
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real vector(3)

        vector = self%BoundFACEs(ID)%moveVector(:)

    end function

    function get_cellCenterOf(self, ID) result(center)
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        integer, intent(in) :: ID
        real center(3)

        center = self%CELLs(ID)%center

    end function

    function get_allOfCellCenters(self) result(centers)
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
        class(FlowFieldUnstructuredGrid), intent(in) :: self
        character(*), intent(in) :: name
        integer info

        select case(name)
        case('node')
            info = size(self%NODEs)
        case('cell')
            info = size(self%CELLs)
        end select

    end function

    
end module unstructuredGrid_m