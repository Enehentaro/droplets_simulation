program CUBE2USG
    use CUBE_mod
    use vtkMesh_operator_m
    use simpleFile_reader
    use caseName_m
    implicit none
    character(100) F_fname, USG_fname, caseName
    character(50),parameter :: filename = 'name.txt'
    character(50), allocatable :: field_name(:)
    character(20), parameter :: CorrespondenceFName = 'vtkCell2cubeNode.bin'
    integer i, j, i_node, n, num_node, num_record, num_cell, nc_max
    integer, pointer :: nc => nowCase
    real X(3)
    real, allocatable :: velocity(:,:)
    logical existance
    type(vtkMesh) USG

    type nodeInfo
        integer cubeID, nodeID(3)
    end type

    type(nodeInfo), allocatable :: vtkCell2cubeNode(:)

    call read_textRecord(filename, field_name)
    num_record = size(field_name)
    
    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call USG%read(USG_fname)

    num_cell = size(USG%cell_array)

    call case_check(num_case = nc_max)
    do nc = 1, nc_max
        caseName = get_caseName()

        do j = 1, num_record
            F_fname = trim(caseName)//'/output/'//field_name(j)

            call read_CUBE_data(F_fname, '')
        
            if (.not.allocated(vtkCell2cubeNode)) then
                inquire(file=CorrespondenceFName, exist=existance)

                if(existance) then
                    call read_nodeInfo

                else
                    allocate(vtkCell2cubeNode(0 : num_cell - 1))

                    do i = 0, num_cell-1
                        X(:) = 0.0
                        num_node = size(USG%cell_array(i)%nodeID)
                        do n = 1, num_node
                            i_node = USG%cell_array(i)%nodeID(n)
                            X(:) = X(:) + USG%node_array(i_node)%coordinate(:)
                        end do
                        X(:) = X(:) / real(num_node)
                        vtkCell2cubeNode(i)%cubeID = get_cube_contains(X)    
                        vtkCell2cubeNode(i)%nodeID(:) = nearest_node(X, vtkCell2cubeNode(i)%cubeID)
                    end do

                    call output_nodeInfo

                end if

            end if

            if (.not.allocated(velocity)) allocate(velocity(3, num_cell))
            do i = 0, num_cell-1
                velocity(:,i+1) = get_velocity_f(vtkCell2cubeNode(i)%nodeID(:), vtkCell2cubeNode(i)%cubeID)
            end do
        
            i = len_trim(F_fname)
            call output_array_asBinary(fname=F_fname(:i-2)//'.array', array=velocity)
            ! call USG%output(F_fname(:i-2)//'.vtk', cellVector=velocity, vectorName='Velocity')
        
        end do

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

        open(newunit=n_unit, file=CorrespondenceFName, form='unformatted', status='old', action='read')
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