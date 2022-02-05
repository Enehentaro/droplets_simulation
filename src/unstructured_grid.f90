module unstructuredGrid_mod
    implicit none
    private

    type node_t
        real coordinate(3)
    end type node_t

    type boundaryTriangle_t
        integer nodeID(3)
        real center(3), normalVector(3), moveVector(3)
    end type boundaryTriangle_t

    type cell_t
        character(5) typeName
        integer, allocatable :: nodeID(:), boundFaceID(:), adjacentCellID(:)
        real center(3), flowVelocity(3), width
    end type cell_t

    type, public :: UnstructuredGrid
        type(node_t), allocatable :: NODEs(:)
        type(cell_t), allocatable :: CELLs(:)

        real MIN_CDN(3), MAX_CDN(3)

        contains

        procedure nearest_cell, nearcell_check, get_MinMaxCDN

        procedure, private :: setupWithFlowFieldFile
        procedure, private :: set_gravity_center, set_MinMaxCDN, point2cellVelocity
        procedure, private :: read_VTK, read_array, read_INP, read_FLD

    end type

    type, public, extends(UnstructuredGrid) :: UnstructuredGridAdjacencySolved
        type(boundaryTriangle_t), allocatable :: BoundFACEs(:)

        integer :: num_refCellSearchFalse = 0, num_refCellSearch = 0

        contains

        procedure updateWithFlowFieldFile
        procedure nearer_cell
        procedure adhesionCheckOnBound
        procedure refCellSearchInfo, search_refCELL

        procedure, private :: AdjacencySolvingProcess
        procedure, private :: read_adjacency, read_boundaries, solve_adacencyOnUnstructuredGrid
        procedure, private :: output_boundaries, output_adjacency, boundary_setting, output_STL

    end type

    public UnstructuredGrid_, UnstructuredGridAdjacencySolved_

    contains

    type(UnstructuredGrid) function UnstructuredGrid_(FlowFieldFile, meshFile)
        character(*), intent(in) :: FlowFieldFile
        character(*), intent(in), optional :: meshFile

        if(present(meshFile)) then
            call UnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile, meshFile)
        else
            call UnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile)
        end if

    end function

    type(UnstructuredGridAdjacencySolved) function UnstructuredGridAdjacencySolved_(FlowFieldFile, meshFile)
        use path_operator_m
        character(*), intent(in) :: FlowFieldFile
        character(*), intent(in), optional :: meshFile
        character(:), allocatable :: Dir

        if(present(meshFile)) then
            UnstructuredGridAdjacencySolved_%UnstructuredGrid = UnstructuredGrid_(FlowFieldFile, meshFile)
        else
            UnstructuredGridAdjacencySolved_%UnstructuredGrid = UnstructuredGrid_(FlowFieldFile)
        end if

        call get_DirFromPath(FlowFieldFile, Dir)

        call UnstructuredGridAdjacencySolved_%AdjacencySolvingProcess(Dir)    !流れ場の前処理

    end function

    subroutine setupWithFlowFieldFile(self, FNAME, meshFile)
        class(UnstructuredGrid) self
        character(*), intent(in) :: FNAME
        character(*), intent(in), optional :: meshFile
        character(:), allocatable :: extension

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME)

            case('inp')
                call self%read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .true., findVelocity = .true.)

                if(.not.allocated(self%CELLs)) then !まだ未割り当てのとき
                    block
                        character(:), allocatable :: topolpgyFNAME
                        topolpgyFNAME = FNAME( : index(FNAME, '_', back=.true.)) // '0' // '.fld' !ゼロ番にアクセス
                        call self%read_FLD(topolpgyFNAME, findTopology= .true., findVelocity = .false.)
                    end block

                    call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)
                end if

            case('array')
                block
                    character(:), allocatable :: FlowDir
                    FlowDir = FNAME( : index(FNAME, '/', back=.true.))
                    call self%read_VTK(FlowDir//meshFile, meshOnly=.true.)
                end block
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', extension
                STOP
                    
        end select

        call self%set_MinMaxCDN()
        call self%set_gravity_center()

    end subroutine

    subroutine updateWithFlowFieldFile(self, FNAME)
        class(UnstructuredGridAdjacencySolved) self
        character(*), intent(in) :: FNAME
        character(:), allocatable :: extension

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME)

            case('inp')
                call self%read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)

            case('array')
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', extension
                STOP
                    
        end select

        call self%set_MinMaxCDN()
        call self%set_gravity_center()

        call self%boundary_setting(first=.false.)
            
    end subroutine

    subroutine AdjacencySolvingProcess(self, dir)
        class(UnstructuredGridAdjacencySolved) self
        character(*), intent(in) :: dir
        logical success

        call self%read_adjacency(dir, success)
        if(success) then
            call self%read_boundaries(dir)

        else
            call self%solve_adacencyOnUnstructuredGrid()
            call self%output_boundaries(dir)
            call self%output_adjacency(dir)

        end if

        call self%boundary_setting(first=.true.)

        call self%output_STL(dir//'shape.stl')

    end subroutine

    subroutine set_MinMaxCDN(self)
        class(UnstructuredGrid) self

        self%MAX_CDN(1) = maxval(self%NODEs(:)%coordinate(1))
        self%MAX_CDN(2) = maxval(self%NODEs(:)%coordinate(2))
        self%MAX_CDN(3) = maxval(self%NODEs(:)%coordinate(3))
        print*, 'MAX_coordinates=', self%MAX_CDN(:)
            
        self%MIN_CDN(1) = minval(self%NODEs(:)%coordinate(1))
        self%MIN_CDN(2) = minval(self%NODEs(:)%coordinate(2))
        self%MIN_CDN(3) = minval(self%NODEs(:)%coordinate(3))
        print*, 'MIN_coordinates=', self%MIN_CDN(:)

    end subroutine

    function get_MinMaxCDN(self) result(MinMax)
        class(UnstructuredGrid) self
        real MinMax(3,2)

        MinMax(:,1) = self%MIN_CDN
        MinMax(:,2) = self%MAX_CDN

    end function

    subroutine read_VTK(self, FNAME, meshOnly)
        use vtkMesh_operator_m
        class(UnstructuredGrid) self
        type(vtkMesh) mesh
        character(*), intent(in) :: FNAME
        logical, intent(in) , optional :: meshOnly
        real, allocatable :: velocity(:,:)
        integer II,KK,IIH, KKMX, IIMX
        logical onlyFlag

        onlyFlag = .false.
        if(present(meshOnly)) then
            if(meshOnly) onlyFlag = .true.
        end if

        if(onlyFlag) then
            call mesh%read(FNAME)
        else
            call mesh%read(FNAME, cellVector=velocity)
        end if

        KKMX = size(mesh%node_array)
        if(.not.allocated(self%NODEs)) allocate(self%NODEs(KKMX))
        do KK = 1, KKMX
            self%NODEs(KK)%coordinate(:) = mesh%node_array(KK-1)%coordinate(:)
        end do
        
        IIMX = size(mesh%cell_array)
        if(.not.allocated(self%CELLs)) allocate(self%CELLs(IIMX))
        do II = 1, IIMX
            IIH = size(mesh%cell_array(II-1)%nodeID)
            self%CELLs(II)%nodeID = mesh%cell_array(II-1)%nodeID(1:IIH) + 1
            select case(mesh%cell_array(II-1)%n_TYPE)
                case(10)
                    self%CELLs(II)%typeName = 'tetra'
                case(13)
                    self%CELLs(II)%typeName = 'prism'
                case(14)
                    self%CELLs(II)%typeName = 'pyrmd'
            end select
        end do

        if(.not.onlyFlag) then
            do II = 1, IIMX
                self%CELLs(II)%flowVelocity(:) = velocity(:,II)
            end do
        endif

        ! print*, NODEs(KKMX)%coordinate(:)
        ! print*, CELLs(IIMX)%flowVelocity(:)
            
    end subroutine

    subroutine read_Array(self, FNAME)
        use simpleFile_reader
        class(UnstructuredGrid) self
        character(*), intent(in) :: FNAME
        real, allocatable :: velocity(:,:)
        integer II

        call read_array_asBinary(FNAME, velocity)

        if(size(self%CELLS) /= size(velocity, dim=2)) then
            print*, 'SIZE ERROR:', size(self%CELLS), size(velocity, dim=2)
            stop
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
        class(UnstructuredGrid) self
        character(*), intent(in) :: FNAME
        INTEGER II,II2,KK,AAmax, IIMX2, AA, n_unit
        integer KKMX, IIMX
        character(6) cellshape
        real, allocatable :: UVWK(:,:)
        integer, allocatable :: ICN2(:,:)
        integer, allocatable :: CELL_TYPE2(:)

        print*, 'READ_INP:', trim(FNAME)
            
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
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
                        self%CELLs(II)%typeName = 'tetra'
                    case(1)
                        self%CELLs(II)%nodeID = ICN2(1:6, II)
                        self%CELLs(II)%typeName = 'prism'
                    case(2)
                        self%CELLs(II)%nodeID = ICN2(1:5, II)
                        self%CELLs(II)%typeName = 'pyrmd'
                end select
            end do
        end if
            
        call self%point2cellVelocity(UVWK)
            
    end subroutine

    subroutine read_FLD(self, FNAME, findTopology, findVelocity)
        use SCT_file_reader_m
        class(UnstructuredGrid) self
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
                    stop
                end if
    
                allocate(self%CELLs(iimx))
                allocate(self%NODEs(kkmx))
    
                do kk = 1, kkmx
                    self%NODEs(kk)%coordinate(:) = real(points(:,kk))
                end do
    
                iicnt = 1
                do ii = 1, iitet
                    self%CELLs(iicnt)%nodeID = tetras(:,ii)
                    self%CELLs(iicnt)%typeName = 'tetra'
                    iicnt = iicnt + 1
                end do
                do ii = 1, iiwed
                    self%CELLs(iicnt)%nodeID = wedges(:,ii)
                    self%CELLs(iicnt)%typeName = 'prism'
                    iicnt = iicnt + 1
                end do
                do ii = 1, iipyr
                    self%CELLs(iicnt)%nodeID = pyramids(:,ii)
                    self%CELLs(iicnt)%typeName = 'pyrmd'
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

    subroutine set_gravity_center(self) !セル重心の算出
        class(UnstructuredGrid) self
        integer II,IIMX, n, num_node, ID

        IIMX = size(self%CELLs)

        !$omp parallel do private(num_node, ID)
        DO II = 1, IIMX
            self%CELLs(II)%center(:) = 0.0
            num_node = size(self%CELLs(II)%nodeID)
            do n = 1, num_node
                ID = self%CELLs(II)%nodeID(n)
                self%CELLs(II)%center(:) = self%CELLs(II)%center(:) + self%NODEs(ID)%coordinate(:)
            end do
            self%CELLs(II)%center(:) = self%CELLs(II)%center(:) / num_node
            
            self%CELLs(II)%width = norm2( &
                self%NODEs(self%CELLs(II)%nodeID(2))%coordinate(:) - self%NODEs(self%CELLs(II)%nodeID(1))%coordinate(:))
        
        END DO
        !$omp end parallel do 
    end subroutine

    subroutine read_adjacency(self, path, success)
        use filename_mod, only : adjacencyFileName
        class(UnstructuredGridAdjacencySolved) self
        character(*), intent(in) :: path
        logical, intent(out) :: success
        integer II,NA, n_unit, num_cells, num_adj, num_BF, NCMAX
        character(:), allocatable :: FNAME
        character(255) str
                
        FNAME = trim(path)//adjacencyFileName
        inquire(file = FNAME, exist = success)
        if(.not.success) return

        print*, 'READ : ', FNAME

        open(newunit=n_unit, FILE=FNAME, STATUS='OLD')
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
        use filename_mod, only : adjacencyFileName
        class(UnstructuredGridAdjacencySolved) self
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
        use filename_mod, only : boundaryFileName
        class(UnstructuredGridAdjacencySolved) self
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(:), allocatable :: FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'READ : ', FNAME
        open(newunit=n_unit, FILE=FNAME , STATUS='old')
            read(n_unit,*) JBMX
            allocate(self%BoundFACEs(JBMX))
            do JB = 1, JBMX
                read(n_unit,*) self%BoundFACEs(JB)%nodeID(:)
            end do
        close(n_unit)
        
    end subroutine

    subroutine output_boundaries(self, path)
        use filename_mod, only : boundaryFileName
        class(UnstructuredGridAdjacencySolved) self
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
        class(UnstructuredGrid) self
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
        
    end function

    integer function nearer_cell(self, X, NCN)  !近セルの探索（隣接セルから）
        class(UnstructuredGridAdjacencySolved) self
        integer, intent(in) :: NCN
        real, intent(in) :: X(3)
        integer NA, featuredCELL, adjacentCELL
        real distance, distance_min
        logical update

        nearer_cell = NCN
        distance_min = norm2(self%CELLs(nearer_cell)%center(:) - X(:))   !注目セル重心と粒子との距離
        update = .true.
        do while(update)    !更新が起こり続ける限り繰り返し
            update = .false.
            featuredCELL = nearer_cell

            checkAdjacent : do NA = 1, size(self%CELLs(featuredCELL)%adjacentCellID)  !全隣接セルに対してループ。

                adjacentCELL = self%CELLs(featuredCELL)%adjacentCellID(NA)  !注目セルの隣接セルのひとつに注目
                if (adjacentCELL <= 0) cycle checkAdjacent

                distance = norm2(self%CELLs(adjacentCELL)%center(:) - X(:))   !注目セル重心と粒子との距離
                if(distance < distance_min) then
                    nearer_cell = adjacentCELL
                    distance_min = distance
                    update = .true.

                end if

            end do checkAdjacent

        end do
        
        ! check:DO
        !     distance(:) = 1.0d10     !初期値はなるべく大きくとる
    
        !     DO NA = 1, size(CELLs(nearer_cell)%adjacentCellID)  !全隣接セルに対してループ。
        !         IIaround = CELLs(nearer_cell)%adjacentCellID(NA)  !現時点で近いとされるセルの隣接セルのひとつに注目
        !         IF (IIaround > 0) distance(NA) = norm2(CENC(:,IIaround)-X(:))   !注目セル重心と粒子との距離を距離配列に代入
        !     END DO
    
        !     distancecheck(2) = minval(distance,dim=1)     !距離配列の最小値
    
        !     if(distancecheck(2) < distancecheck(1)) then !より近いセルの発見で条件満足
        !         distancecheck(1) = distancecheck(2)    !最小値の更新
        !         index_min = minloc(distance,dim=1)            !最小値のインデックス
        !         nearer_cell = NEXT_CELL(index_min, nearer_cell)    !現時点で近いとされるセルの更新
        !         if(nearer_cell==0) then
        !             print*,'nearer_cell_error', nearer_cell, X(:)
        !             return
        !         end if

        !     else  !より近いセルを発見できなかった場合
    
        !         exit check     !ループ脱出
    
        !     end if
        
        ! END DO check
        
    end function

    logical function nearcell_check(self, X, NCN)
        class(UnstructuredGrid) self
        real, intent(in) :: X(3)
        integer, intent(in) :: NCN
        real :: distance

        distance = norm2(X(:) - self%CELLs(NCN)%center(:))

        if (distance < 1.0d1*self%CELLs(NCN)%width) then
            nearcell_check = .True.
        else
            nearcell_check = .False.
            ! print*, 'nearcell_check:False', distance, CELLs(NCN)%width
        end if

    end function

    subroutine search_refCELL(self, X, reference_cell, stat)
        class(UnstructuredGridAdjacencySolved) self
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
        class(UnstructuredGridAdjacencySolved) self
        character(*), intent(in) :: name

        select case(name)
            case('NumFalse')
                refCellSearchInfo = self%num_refCellSearchFalse
            case('FalseRate')
                refCellSearchInfo = 100 * self%num_refCellSearchFalse / (self%num_refCellSearch + 1)
            case default
                print*, 'ERROR refCellSearchInfo : ', name
                stop
        end select

    end function
                     
    subroutine boundary_setting(self, first) !全境界面に対して外向き法線ベクトルと重心を算出
        use vector_m
        class(UnstructuredGridAdjacencySolved) self
        logical, intent(in) :: first
        integer II, JJ, JB, IIMX, JBMX, nodeID(3)
        real :: a(3), b(3), r(3), normalVector(3)
        type(boundaryTriangle_t), allocatable :: BoundFACEs_pre(:)

        IIMX = size(self%CELLs)

        if(.not.first) BoundFACEs_pre = self%BoundFACEs
        
        print*, 'SET:boundary'
        
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
        class(UnstructuredGridAdjacencySolved) self
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
        class(UnstructuredGrid) self
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
        class(UnstructuredGridAdjacencySolved) self
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

    subroutine solve_adacencyOnUnstructuredGrid(self)
        use adjacencySolver_m
        class(UnstructuredGridAdjacencySolved) self
        integer i, j, num_adjacent, num_boundFace
        integer, parameter :: max_vertex=6, max_adjacent=4, max_boundFace=4
        integer cellVertices(max_vertex, size(self%CELLs))
        integer adjacentCellArray(max_adjacent, size(self%CELLs))
        integer cellBoundFaces(max_boundFace, size(self%CELLs))
        integer, allocatable :: boundFaceVertices(:,:)

        cellVertices = 0
        do i = 1, size(self%CELLs)
            cellVertices(1:size(self%CELLs(i)%nodeID(:)), i) = self%CELLs(i)%nodeID(:)
        end do

        cellBoundFaces = 0
        adjacentCellArray = 0
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

    integer function get_mesh_info(self, name)
        class(UnstructuredGrid) self
        character(*), intent(in) :: name
        
        select case(name)
            case('node')
                get_mesh_info = size(self%NODEs)

            case('cell')
                get_mesh_info = size(self%CELLs)

            case('tetra')
                get_mesh_info = count(self%CELLs(:)%typeName == 'tetra')

            case('prism')
                get_mesh_info = count(self%CELLs(:)%typeName == 'prism')

            case('pyramid')
                get_mesh_info = count(self%CELLs(:)%typeName == 'pyrmd')

            case default
                get_mesh_info = -1

        end select

    end function

    function extensionOf(FileName) result(extension)
        character(*), intent(in) :: FileName
        character(:), allocatable :: extension

        extension = FileName(index(FileName, '.')+1 : )

    end function
      
    ! subroutine deallocation_unstructuredGRID
    !     use vtkMesh_operator_m

    !     print*,  '**Deallocation** : unstructuredGRID'

    !     deallocate(CELLs)
    !     deallocate(NODEs)
    !     deallocate(BoundFACEs)

    ! end subroutine
    
end module unstructuredGrid_mod