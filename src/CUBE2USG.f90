program CUBE2USG
    use CUBE_mod
    use vtkMesh_operator_m
	use counter
    implicit none
    character(50) F_fname, USG_fname
    character(50), allocatable :: field_name(:)
    integer i, j, n, num_node, num_file, num_cell
    real X(3)

    type nodeInfo
        integer cubeID, nodeID(3)
    end type

    type(nodeInfo), allocatable :: vtkCell2cubeNode(:)

    num_file=count_txt()
    print *, num_file

    allocate(field_name(num_file))
    call read_txt(num_file,field_name)
    
    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call read_VTK_mesh(USG_fname, meshONLY=.true.)

    num_cell = size(cell_array)

    do j = 1, num_file
        F_fname = field_name(j)

        call read_CUBE_data(F_fname, '')
    
        if (.not.allocated(vtkCell2cubeNode)) then
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

        end if

        do i = 0, num_cell-1
            cell_array(i)%vector(:) = get_velocity_f(vtkCell2cubeNode(i)%nodeID(:), vtkCell2cubeNode(i)%cubeID)
        end do
    
        i = len_trim(F_fname)
        call output_VTK_mesh(F_fname(:i-2)//'.vtk')
    
    end do
    
end program CUBE2USG