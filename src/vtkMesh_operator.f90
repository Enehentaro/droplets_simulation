module vtkMesh_operator_m
    implicit none
    private

    type node_onVTK_t
        real coordinate(3)
    end type

    type cell_onVTK_t
        integer, allocatable :: nodeID(:)
        integer n_TYPE
    end type

    type, public :: vtkMesh
        type(node_onVTK_t), private, allocatable :: node_array(:)
        type(cell_onVTK_t), private, allocatable :: cell_array(:)

        contains

        procedure, private :: allocation_node, allocation_cell
        procedure :: read => read_vtkMesh
        procedure :: output => output_vtkMesh
        procedure get_numCell, get_numNode, get_nodeCoordinate, get_cellVertices, get_cellCenter
        procedure set_nodeCoordinate, set_cellVertices

    end type

    public vtkMesh_

    contains

    type(vtkMesh) function vtkMesh_(cdn, vertices, types)
        real, intent(in) :: cdn(:,:)
        integer, intent(in) :: vertices(:,:), types(:)

        call vtkMesh_%set_nodeCoordinate(cdn)
        call vtkMesh_%set_cellVertices(vertices, types)

    end function

    subroutine read_vtkMesh(self, FNAME, action, cellScalar, cellVector)
        class(vtkMesh) self
        character(*), intent(in) :: FNAME
        character(*), intent(in), optional ::  action
        real, allocatable, intent(out), optional :: cellScalar(:), cellVector(:,:)
        integer II,KK,IIH, n_unit, KKMX, IIMX, ios
        character AAA*7, str*99
        logical dataOnly

        dataOnly = .false.
        if(present(action)) then
            if(action=='dataOnly') dataOnly = .true.
        end if

        print*, 'READ_VTK:', FNAME
            
        open(newunit=n_unit,FILE=FNAME, status='old', action='read')
            if(dataOnly) then
                do while(index(AAA, 'CELL_DATA')==0)
                    read(n_unit, '(A)') AAA
                end do
            else
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,*) AAA,KKMX
                    
                call self%allocation_node(KKMX)
                DO KK = 0, KKMX-1
                    read(n_unit,*) self%node_array(KK)%coordinate(:)
                END DO
                read(n_unit,'()')
                read(n_unit, *) AAA,IIMX
                
                call self%allocation_cell(IIMX)
       
                DO II = 0, IIMX-1
                    read(n_unit, '(A)') str !一旦1行丸ごと読む
                    read(str, *) IIH
                    allocate(self%cell_array(II)%nodeID(IIH))
                    read(str, *) IIH, self%cell_array(II)%nodeID(:)
                    self%cell_array(II)%nodeID(:) = self%cell_array(II)%nodeID(:)
                END DO

                read(n_unit,'()')

                read(n_unit,'()')  !CELL_TYPES
                DO II = 0, IIMX-1
                    read(n_unit, *) self%cell_array(II)%n_TYPE
                END DO
                if((.not.present(cellScalar)) .and. (.not.present(cellVector))) return
                read(n_unit,'()')
    
                read(n_unit,'()')   !CELL_DATA

            end if

            do 
                read(n_unit,'(A)', iostat=ios) AAA
                if(ios/=0) exit
                select case(AAA)
                    case('SCALARS')
                        read(n_unit,'()')
                        if(present(cellScalar)) then
                            allocate(cellScalar(IIMX))
                            DO II = 1, IIMX
                                read(n_unit, *) cellScalar(II)
                            END DO
                        else
                            DO II = 1, IIMX
                                read(n_unit, '()')
                            END DO
                        end if
                    case('VECTORS')
                        if(present(cellVector)) then
                            allocate(cellVector(3,IIMX))
                            DO II = 1, IIMX
                                read(n_unit, *) cellVector(:,II)
                            END DO
                        end if
                        exit
                end select
            end do
            
        close(n_unit)
            
    end subroutine

    subroutine output_vtkMesh(self, FNAME, cellScalar, cellVector, scalarName, vectorName)
        class(vtkMesh) self
        character(*), intent(in) :: FNAME
        real, intent(in), optional :: cellScalar(:), cellVector(:,:)
        character(*), intent(in), optional :: scalarName, vectorName
        integer II,KK, n_unit, KKMX, IIMX, IITOTAL

        print*, 'OUTPUT_VTK:', FNAME
            
        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit, '(A)') '# vtk DataFile Version 2.0'
            write(n_unit, '(A)') 'Header'
            write(n_unit, '(A)') 'ASCII'
            write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

            KKMX = size(self%node_array)
            write(n_unit, '(A,1x,I0,1x,A)') 'POINTS', KKMX, 'float'   
            DO KK = 0, KKMX-1
                write(n_unit,'(3(e12.5,2X))') self%node_array(KK)%coordinate(:)
            END DO
            write(n_unit,'()')

            IIMX = size(self%cell_array)
            IITOTAL = 0
            do II = 0, IIMX-1
                IITOTAL = IITOTAL + size(self%cell_array(II)%nodeID)  + 1
            end do
            write(n_unit,'(A,I0,2X,I0)') 'CELLS ', IIMX, IITOTAL
            DO II = 0, IIMX-1
                select case(self%cell_array(II)%n_TYPE)
                    case(10)
                        write(n_unit, '(5(I0,2X))') 4, self%cell_array(II)%nodeID(1:4)
                    case(11)
                        write(n_unit, '(9(I0,2X))') 8, self%cell_array(II)%nodeID(1:8)
                    case(13)
                        write(n_unit, '(7(I0,2X))') 6, self%cell_array(II)%nodeID(1:6)
                    case(14)
                        write(n_unit, '(6(I0,2X))') 5, self%cell_array(II)%nodeID(1:5)
                end select
            END DO
            write(n_unit,'()')

            write(n_unit,'(A,I0)') 'CELL_TYPES ', IIMX
            DO II = 0, IIMX-1
                write(n_unit, '(I0)') self%cell_array(II)%n_TYPE
            END DO

            if(present(cellScalar) .or. present(cellVector)) then
                write(n_unit,'()')
                write(n_unit,'(A,I0)') 'CELL_DATA ', IIMX

                if(present(cellScalar)) then
                    if(present(scalarName)) then
                        write(n_unit,'(A)') 'SCALARS '//scalarName//' float'
                    else
                        write(n_unit,'(A)') 'SCALARS scalar float'
                    end if
                    write(n_unit,'(A)') 'LOOKUP_TABLE default'
                    DO II = 1, IIMX
                        write(n_unit,'(e12.5)') cellScalar(II)
                    END DO
                end if

                if(present(cellVector)) then
                    if(present(vectorName)) then
                        write(n_unit,'(A)') 'VECTORS '//vectorName//' float'
                    else
                        write(n_unit,'(A)') 'VECTORS vector float' 
                    end if  
                    DO II = 1, IIMX
                        write(n_unit,'(3(e12.5,2X))') cellVector(:, II)
                    END DO
                end if

            end if
            
        close(n_unit)
            
    end subroutine

    subroutine allocation_node(self, num_node)
        class(vtkMesh) self
        integer, intent(in) :: num_node
        if(allocated(self%node_array))  deallocate(self%node_array)
        allocate(self%node_array(0 : num_node-1))     
    end subroutine

    subroutine allocation_cell(self, num_cell)
        class(vtkMesh) self
        integer, intent(in) :: num_cell
        if(allocated(self%cell_array))  deallocate(self%cell_array)
        allocate(self%cell_array(0 : num_cell-1))        
    end subroutine

    integer function get_numNode(self)
        class(vtkMesh), intent(in) :: self
        get_numNode = size(self%node_array)
    end function

    integer function get_numCell(self)
        class(vtkMesh), intent(in) :: self
        get_numCell = size(self%cell_array)
    end function

    function get_nodeCoordinate(self) result(cdn)
        class(vtkMesh), intent(in) :: self
        real, allocatable :: cdn(:,:)
        integer k, num_node

        num_node = size(self%node_array)
        allocate(cdn(3, num_node))

        do k = 1, num_node
            cdn(:,k) = self%node_array(k-1)%coordinate(:)
        end do

    end function

    subroutine get_cellVertices(self, vertices, types)
        class(vtkMesh), intent(in) :: self
        integer, allocatable, intent(out) :: vertices(:,:), types(:)
        integer i, num_cell, num_node

        num_cell = size(self%cell_array)
        allocate(vertices(6, num_cell))
        allocate(types(num_cell))

        do i = 1, num_cell
            num_node = size(self%cell_array(i-1)%nodeID)
            vertices(1:num_node, i) = self%cell_array(i-1)%nodeID(:) + 1
            types(i) = self%cell_array(i-1)%n_TYPE
        end do

    end subroutine

    subroutine set_nodeCoordinate(self, cdn)
        class(vtkMesh) self
        real, intent(in) :: cdn(:,:)
        integer k, num_node

        num_node = size(cdn, dim=2)
        call self%allocation_node(num_node)

        do k = 1, num_node
            self%node_array(k-1)%coordinate(:) = cdn(:,k) 
        end do

    end subroutine

    subroutine set_cellVertices(self, vertices, types)
        class(vtkMesh) self
        integer, intent(in) :: vertices(:,:), types(:)
        integer i, num_cell, num_node

        num_cell = size(vertices, dim=2)
        call self%allocation_cell(num_cell)

        do i = 1, num_cell
            select case(types(i))
                case(10)!tetra
                    num_node = 4
                case(11)!hexa
                    num_node = 8
                case(13)!prism
                    num_node = 6
                case(14)!pyramid
                    num_node = 5
            end select
            self%cell_array(i-1)%nodeID = vertices(1:num_node, i) - 1
            self%cell_array(i-1)%n_TYPE = types(i)
        end do

    end subroutine

    function get_cellCenter(self) result(center)
        class(vtkMesh), intent(in) :: self
        real, allocatable :: center(:,:)
        real x(3)
        integer i, k, num_node, num_cell, nodeID

        num_cell = size(self%cell_array)
        allocate(center(3, num_cell))

        do i = 1, num_cell
            x(:) = 0.0
            num_node = size(self%cell_array(i-1)%nodeID)
            do k = 1, num_node
                nodeID = self%cell_array(i-1)%nodeID(k)
                x(:) = x(:) + self%node_array(nodeID)%coordinate(:)
            end do
            center(:,) = x(:) / real(num_node)
        end do

    end function
    
end module vtkMesh_operator_m
