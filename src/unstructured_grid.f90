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
        real center(3), flowVelocity(3), threshold
    end type cell_t

    type, public :: UnstructuredGrid
        type(node_t), allocatable :: NODEs(:)
        type(cell_t), allocatable :: CELLs(:)
        type(boundaryTriangle_t), allocatable :: BoundFACEs(:)

        real, private :: MIN_CDN(3), MAX_CDN(3)
        integer, private :: num_refCellSearchFalse = 0, num_refCellSearch = 0

        contains

        procedure nearest_cell, nearcell_check, get_MinMaxCDN

        procedure, private :: setupWithFlowFieldFile
        procedure, private :: set_gravity_center, set_MinMaxCDN, point2cellVelocity
        procedure, private :: read_VTK, read_array, read_INP, read_FLD

        !=====================================================================

        procedure updateWithFlowFieldFile
        procedure nearer_cell
        procedure adhesionCheckOnBound
        procedure refCellSearchInfo, search_refCELL

        procedure, private :: AdjacencySolvingProcess
        procedure, private :: read_adjacency, read_boundaries, solve_adacencyOnUnstructuredGrid
        procedure, private :: output_boundaries, output_adjacency, boundary_setting, output_STL

    end type

    public UnstructuredGrid_

    contains

    type(UnstructuredGrid) function UnstructuredGrid_(FlowFieldFile, meshFile)
        use path_operator_m
        character(*), intent(in) :: FlowFieldFile
        character(*), intent(in), optional :: meshFile
        character(:), allocatable :: Dir

        if(present(meshFile)) then
            call UnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile, meshFile)
            call get_DirFromPath(meshFile, Dir)
        else
            call UnstructuredGrid_%setupWithFlowFieldFile(FlowFieldFile)
            call get_DirFromPath(FlowFieldFile, Dir)
        end if

        call UnstructuredGrid_%AdjacencySolvingProcess(Dir)    !?????????????????????

    end function

    subroutine setupWithFlowFieldFile(self, FNAME, meshFile)
        class(UnstructuredGrid) self
        character(*), intent(in) :: FNAME
        character(*), intent(in), optional :: meshFile
        character(:), allocatable :: extension

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME, meshOnly=.false.)

            case('inp')
                call self%read_INP(FNAME)   !INP???????????????(SHARP???)

            case('fld')
                call self%read_FLD(FNAME, findTopology= .true., findVelocity = .true.)

                if(.not.allocated(self%CELLs)) then !??????????????????????????????
                    block
                        character(:), allocatable :: topolpgyFNAME
                        topolpgyFNAME = FNAME( : index(FNAME, '_', back=.true.)) // '0' // '.fld' !????????????????????????
                        call self%read_FLD(topolpgyFNAME, findTopology= .true., findVelocity = .false.)
                    end block

                    call self%read_FLD(FNAME, findTopology= .false., findVelocity = .true.)
                end if

            case('array')
                call self%read_VTK(meshFile, meshOnly=.true.)
                call self%read_Array(FNAME)

            case default
                print*,'FILE_EXTENSION NG : ', extension
                STOP
                    
        end select

        call self%set_MinMaxCDN()
        call self%set_gravity_center()

    end subroutine

    subroutine updateWithFlowFieldFile(self, FNAME)
        class(UnstructuredGrid) self
        character(*), intent(in) :: FNAME
        character(:), allocatable :: extension

        select case(extensionOf(FNAME))
            case('vtk')
                call self%read_VTK(FNAME, meshOnly=.false.)

            case('inp')
                call self%read_INP(FNAME)   !INP???????????????(SHARP???)

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
        class(UnstructuredGrid) self
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
        logical, intent(in) :: meshOnly
        real, allocatable :: cdn(:,:)
        integer, allocatable :: vertices(:,:), types(:)
        real, allocatable :: velocity(:,:)
        integer II,KK, KKMX, IIMX

        if(meshOnly) then
            call mesh%read(FNAME)
        else
            call mesh%read(FNAME, cellVector=velocity)
        end if

        KKMX = mesh%get_numNode()
        if(.not.allocated(self%NODEs)) allocate(self%NODEs(KKMX))
        cdn = mesh%get_nodeCoordinate()
        do KK = 1, KKMX
            self%NODEs(KK)%coordinate(:) = cdn(:,KK)
        end do
        
        IIMX = mesh%get_numCell()
        if(.not.allocated(self%CELLs)) allocate(self%CELLs(IIMX))
        call mesh%get_cellVertices(vertices, types)
        do II = 1, IIMX
            select case(types(II))
                case(10)
                    self%CELLs(II)%nodeID = vertices(1:4, II)
                    self%CELLs(II)%typeName = 'tetra'
                case(13)
                    self%CELLs(II)%nodeID = vertices(1:6, II)
                    self%CELLs(II)%typeName = 'prism'
                case(14)
                    self%CELLs(II)%nodeID = vertices(1:5, II)
                    self%CELLs(II)%typeName = 'pyrmd'
            end select
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
        !  INP??????????????????????????????????????????????????????????????????????????????
        class(UnstructuredGrid) self
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
    
                read(n_unit,fmt='(I10)',advance='no')AA  !????????????????????????????????????
                read(n_unit,fmt='(I6)',advance='no')AA  !???????????????????????????????????????????????????
                read(n_unit,fmt='(A6)',advance='no')cellshape
    
                cellshape = adjustl(cellshape)    !?????????
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
                        read(n_unit,*)ICN2(5,II),ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II) !INP?????????????????????????????????VTK?????????????????????????????????????????????????????????????????????
                        if (ICN2(1,II)==0.or.ICN2(5,II)==0) print*, 'ICN2_WARNING_pyr:', ICN2(:,II)
        
                    end if
    
                ELSE
                    read(n_unit,'()')  !???????????????????????????????????????????????????????????????????????????
    
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
            !!????????????????????????, ??????????????????????????????????????????????????????????????????????????????.  
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

    subroutine set_gravity_center(self) !??????????????????????????????????????????
        class(UnstructuredGrid) self
        integer II,IIMX, n, num_node, nodeID
        real vector(3)

        IIMX = size(self%CELLs)

        !$omp parallel do private(num_node, nodeID, vector)
        DO II = 1, IIMX
            vector(:) = 0.0
            num_node = size(self%CELLs(II)%nodeID)
            do n = 1, num_node
                nodeID = self%CELLs(II)%nodeID(n)
                vector(:) = vector(:) + self%NODEs(nodeID)%coordinate(:)
            end do
            self%CELLs(II)%center(:) = vector(:) / real(num_node)
            
            self%CELLs(II)%threshold = 0.0
            do n = 1, num_node
                nodeID = self%CELLs(II)%nodeID(n)
                vector(:) = self%NODEs(nodeID)%coordinate(:) - self%CELLs(II)%center(:)
                self%CELLs(II)%threshold = max(self%CELLs(II)%threshold, norm2(vector)) 
            end do
            ! print*, self%CELLs(II)%threshold
        
        END DO
        !$omp end parallel do

    end subroutine

    subroutine read_adjacency(self, path, success)
        use filename_mod, only : adjacencyFileName
        class(UnstructuredGrid) self
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
        use filename_mod, only : adjacencyFileName
        class(UnstructuredGrid) self
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
        class(UnstructuredGrid) self
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
        use filename_mod, only : boundaryFileName
        class(UnstructuredGrid) self
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

    integer function nearest_cell(self, X) !??????????????????NCN?????????
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
        
        nearest_cell = minloc(distance, dim=1)   !???????????????????????????
        ! print*, distance(nearest_cell)
        
    end function

    integer function nearer_cell(self, X, NCN)  !??????????????????????????????????????????
        class(UnstructuredGrid) self
        integer, intent(in) :: NCN
        real, intent(in) :: X(3)
        integer NA, featuringCellID, adjaCellID
        integer, allocatable :: adjacentCellIDs(:)
        real distance, distance_min
        logical update

        nearer_cell = NCN
        distance_min = norm2(self%CELLs(nearer_cell)%center(:) - X(:))   !???????????????????????????????????????
        update = .true.

        do while(update)    !?????????????????????????????????????????????
            update = .false.

            featuringCellID = nearer_cell       !????????????

            adjacentCellIDs = self%CELLs(featuringCellID)%adjacentCellID(:)  !??????????????????????????????

            do NA = 1, size(adjacentCellIDs)  !??????????????????????????????????????????????????????

                adjaCellID = adjacentCellIDs(NA)  !?????????????????????????????????????????????
                ! if (adjacentCELL <= 0) cycle checkAdjacent

                distance = norm2(self%CELLs(adjaCellID)%center(:) - X(:))   !???????????????????????????????????????
                if(distance < distance_min) then
                    nearer_cell = adjaCellID
                    distance_min = distance
                    update = .true.
                end if

            end do

        end do
        
    end function

    logical function nearcell_check(self, X, NCN)
        class(UnstructuredGrid) self
        real, intent(in) :: X(3)
        integer, intent(in) :: NCN
        real :: distance

        distance = norm2(X(:) - self%CELLs(NCN)%center(:))

        !??????????????????????????????????????????????????????????????????
        !?????????????????????????????????????????????????????????OK?????????????????????????????????????????????
        if (distance < self%CELLs(NCN)%threshold) then
            nearcell_check = .True.
        else
            nearcell_check = .False.
            ! print*, 'nearcell_check:False', distance, self%CELLs(NCN)%threshold
        end if

    end function

    subroutine search_refCELL(self, X, reference_cell, stat)
        class(UnstructuredGrid) self
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
        class(UnstructuredGrid) self
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
                     
    subroutine boundary_setting(self, first) !?????????????????????????????????????????????????????????????????????
        use vector_m
        class(UnstructuredGrid) self
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
            
                r(:) = self%CELLs(II)%center(:) - self%BoundFACEs(JB)%center(:)  !?????????????????????????????????????????????
                if(dot_product(normalVector(:), r(:)) > 0.0) then
                    normalVector(:) = normalVector(:) * (-1.0) !??????????????????????????????????????????????????????
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
        class(UnstructuredGrid) self
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
            !??????????????????????????????????????????????????????????????????????????????????????????????????????

            if (inner + radius > 0.d0) then !(???????????????+????????????)?????????????????????????????????  
                stat = JB       !???????????????????????????
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
        class(UnstructuredGrid) self
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
        class(UnstructuredGrid) self
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

        extension = FileName(index(FileName, '.', back=.true.)+1 : )

    end function
      
    ! subroutine deallocation_unstructuredGRID
    !     use vtkMesh_operator_m

    !     print*,  '**Deallocation** : unstructuredGRID'

    !     deallocate(CELLs)
    !     deallocate(NODEs)
    !     deallocate(BoundFACEs)

    ! end subroutine
    
end module unstructuredGrid_mod