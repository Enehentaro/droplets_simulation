module vtkMesh_operator_m
    implicit none

    type, private :: node_onVTK_t
        real coordinate(3)
        real scalar
        real vector(3)
    end type node_onVTK_t

    type, private :: cell_onVTK_t
        integer, allocatable :: nodeID(:)
        integer n_TYPE
        real scalar
        real vector(3)
    end type cell_onVTK_t

    type(node_onVTK_t), allocatable :: node_array(:)
    type(cell_onVTK_t), allocatable :: cell_array(:)

    contains

    subroutine read_VTK_mesh(FNAME, meshONLY)
        character(*), intent(in) :: FNAME
        logical, optional ::  meshONLY
        integer II,KK,IIH, n_unit, KKMX, IIMX, ios
        character AAA*7, str*99

        print*, 'READ_VTK:', FNAME
            
        open(newunit=n_unit,FILE=FNAME, STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,*) AAA,KKMX
                
            call allocation_node(KKMX)
            DO KK = 0, KKMX-1
                read(n_unit,*) node_array(KK)%coordinate(:)
            END DO
            read(n_unit,'()')
            read(n_unit, *) AAA,IIMX
            
            if(.not.allocated(cell_array)) then
                call allocation_cell(IIMX)
   
                DO II = 0, IIMX-1
                    read(n_unit, '(A)') str !一旦1行丸ごと読む
                    read(str, *) IIH
                    allocate(cell_array(II)%nodeID(IIH))
                    read(str, *) IIH, cell_array(II)%nodeID(:)
                    cell_array(II)%nodeID(:) = cell_array(II)%nodeID(:)
                END DO

            else
                DO II = 0, IIMX-1
                    read(n_unit, '()')
                END DO
                
            end if

            read(n_unit,'()')

            read(n_unit,'()')  !CELL_TYPES
            DO II = 0, IIMX-1
                read(n_unit, *) cell_array(II)%n_TYPE
            END DO
            read(n_unit,'()')

            if(present(meshONLY)) then
                if(meshONLY) return
            end if

            read(n_unit,'()')   !CELL_DATA
            do 
                read(n_unit,'(A)', iostat=ios) AAA
                if(ios/=0) exit
                select case(AAA)
                    case('SCALARS')
                        read(n_unit,'()')
                        DO II = 0, IIMX-1
                            read(n_unit, *) cell_array(II)%scalar
                        END DO
                    case('VECTORS')
                        DO II = 0, IIMX-1
                            read(n_unit, *) cell_array(II)%vector(:)
                        END DO
                        exit
                end select
            end do
            
        close(n_unit)
            
    end subroutine read_VTK_mesh

    subroutine output_VTK_mesh(FNAME, meshONLY)
        character(*), intent(in) :: FNAME
        logical, optional, intent(in) :: meshONLY
        integer II,KK, n_unit, KKMX, IIMX, IITOTAL

        print*, 'OUTPUT_VTK:', FNAME
            
        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit, '(A)') '# vtk DataFile Version 2.0'
            write(n_unit, '(A)') 'Header'
            write(n_unit, '(A)') 'ASCII'
            write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

            KKMX = size(node_array)
            write(n_unit, '(A,1x,I0,1x,A)') 'POINTS', KKMX, 'float'   
            DO KK = 0, KKMX-1
                write(n_unit,'(3(e12.5,2X))') node_array(KK)%coordinate(:)
            END DO
            write(n_unit,'()')

            IIMX = size(cell_array)
            IITOTAL = 5*count(cell_array(:)%n_TYPE==10) + 7*count(cell_array(:)%n_TYPE==13) &
                    + 6*count(cell_array(:)%n_TYPE==14) + 9*count(cell_array(:)%n_TYPE==11)
            write(n_unit,'(A,I0,2X,I0)') 'CELLS ', IIMX, IITOTAL
            DO II = 0, IIMX-1
                select case(cell_array(II)%n_TYPE)
                    case(10)
                        write(n_unit, '(5(I0,2X))') 4, cell_array(II)%nodeID(1:4)
                    case(11)
                        write(n_unit, '(9(I0,2X))') 8, cell_array(II)%nodeID(1:8)
                    case(13)
                        write(n_unit, '(7(I0,2X))') 6, cell_array(II)%nodeID(1:6)
                    case(14)
                        write(n_unit, '(6(I0,2X))') 5, cell_array(II)%nodeID(1:5)
                end select
            END DO
            write(n_unit,'()')

            write(n_unit,'(A,I0)') 'CELL_TYPES ', IIMX
            DO II = 0, IIMX-1
                write(n_unit, '(I0)') cell_array(II)%n_TYPE
            END DO
            if(present(meshONLY)) then
                if(meshONLY) return
            end if
            write(n_unit,'()')

            write(n_unit,'(A,I0)') 'CELL_DATA ', IIMX  
            write(n_unit,'(A)') 'VECTORS Velocity float'    
            DO II = 0, IIMX-1
                write(n_unit,'(3(e12.5,2X))') cell_array(II)%vector(:)
            END DO
            
        close(n_unit)
            
    end subroutine output_VTK_mesh

    subroutine allocation_node(num_node)
        integer, intent(in) :: num_node
        if(.not.allocated(node_array)) allocate(node_array(0 : num_node-1))     
    end subroutine

    subroutine allocation_cell(num_cell)
        integer, intent(in) :: num_cell
        if(.not.allocated(cell_array)) allocate(cell_array(0 : num_cell-1))     
    end subroutine
      
    subroutine deallocation_VTK
        deallocate(node_array)
        deallocate(cell_array)
    end subroutine deallocation_VTK
    
end module vtkMesh_operator_m