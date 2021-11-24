program CUBE2USG
    use CUBE_mod
    use vtkMesh_operator_m
    implicit none
    character(50) F_fname, USG_fname
    character(50), allocatable :: fname(:)
    integer i, j, n, num_node, num_file, num_cell
    real X(3)

    type nodeInfo
        integer cubeID, nodeID(3)
    end type

    type(nodeInfo), allocatable :: nodeInfos(:)

    print *, 'num_file ?'
    read(5,*) num_file

    allocate(fname(num_file))
    
    print *, 'F_FileName ?'
    do i = 1, num_file
        read(5,*) fname(i)
    end do
    
    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call read_VTK_mesh(USG_fname, meshONLY=.true.)

    num_cell = size(cell_array)

    do j = 1, num_file
        F_fname = fname(j)

        call read_CUBE_data(F_fname, '')
    
        if (.not.allocated(nodeInfos)) then
            allocate(nodeInfos(0 : num_cell - 1))

            do i = 0, num_cell-1
                X(:) = 0.0
                num_node = size(cell_array(i)%nodeID)
                do n = 1, num_node
                    X(:) = X(:) + node_array(cell_array(i)%nodeID(n))%coordinate(:)
                end do
                X(:) = X(:) / real(num_node)
                nodeInfos(i)%cubeID = get_cube_contains(X)    
                nodeInfos(i)%nodeID(:) = nearest_node(X, nodeInfos(i)%cubeID)
            end do

        end if

        do i = 0, num_cell-1
            cell_array(i)%vector(:) = get_velocity_f(nodeInfos(i)%nodeID(:), nodeInfos(i)%cubeID)
        end do
    
        i = len_trim(F_fname)
        call output_VTK_mesh(F_fname(:i-2)//'.vtk')
    
    end do
    
end program CUBE2USG