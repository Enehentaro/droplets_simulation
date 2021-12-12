program CUBE2USG
    use CUBE_mod
    use vtkMesh_operator_m
    use csv_reader
    implicit none
    character(50) F_fname, USG_fname
    character(50), allocatable :: SG_fname(:,:)
    character(20), parameter :: CorrespondenceFName = 'vtkCell2cubeNode.bin'
    integer i, j, n, num_node, num_cell
    integer :: num_file(2)
    real X(3)
    logical existance

    type nodeInfo
        integer cubeID, nodeID(3)
    end type

    type(nodeInfo), allocatable :: vtkCell2cubeNode(:)


    call read_csv_char('name.txt', SG_fname, header=.false., mat_size=num_file)

    
    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call read_VTK_mesh(USG_fname, meshONLY=.true.)

    num_cell = size(cell_array)

    do j = 1, num_file(2)
        F_fname = SG_fname(num_file(1),j)

        call read_CUBE_data(F_fname, '')
    
        if (.not.allocated(vtkCell2cubeNode)) then
            inquire(file=CorrespondenceFName, exist=existance)

            if(existance) then
                call read_nodeInfo

            else
                allocate(vtkCell2cubeNode(0 : num_cell - 1))

                do i = 0, num_cell-1
                    X(:) = 0.0
                    num_node = size(cell_array(i)%nodeID)
                    do n = 1, num_node
                        X(:) = X(:) + node_array(cell_array(i)%nodeID(n))%coordinate(:)
                    end do
                    X(:) = X(:) / real(num_node)
                    vtkCell2cubeNode(i)%cubeID = get_cube_contains(X)    
                    vtkCell2cubeNode(i)%nodeID(:) = nearest_node(X, vtkCell2cubeNode(i)%cubeID)
                end do

                call output_nodeInfo

            end if

        end if

        do i = 0, num_cell-1
            cell_array(i)%vector(:) = get_velocity_f(vtkCell2cubeNode(i)%nodeID(:), vtkCell2cubeNode(i)%cubeID)
        end do
    
        i = len_trim(F_fname)
        call output_VTK_mesh(F_fname(:i-2)//'.vtk')
    
    end do

    contains

    subroutine output_nodeInfo
        integer n_unit, k

        print*, 'output: ', CorrespondenceFName

        open(newunit=n_unit, file=CorrespondenceFName, form='unformatted', status='new')
            write(n_unit) num_cell, get_numCube()

            do k = 0, num_cell-1
                write(n_unit) vtkCell2cubeNode(k)
            end do

        close(n_unit)

    end subroutine

    subroutine read_nodeInfo
        integer n_unit, k, num_cell_, num_cube

        print*, 'read: ', CorrespondenceFName

        open(newunit=n_unit, file=CorrespondenceFName, form='unformatted', status='old')
            read(n_unit) num_cell_, num_cube

            if((num_cell_/=num_cell).or.(num_cube/=get_numCube())) then
                print*, 'SizeERROR:', num_cell_, num_cell, num_cube, get_numCube()
                stop
            end if

            allocate(vtkCell2cubeNode(0 : num_cell_ - 1))

            do k = 0, num_cell_-1
                read(n_unit) vtkCell2cubeNode(k)
            end do

        close(n_unit)

    end subroutine
    
end program CUBE2USG