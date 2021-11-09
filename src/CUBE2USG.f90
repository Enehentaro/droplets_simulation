include 'vtkMesh_operator.f90'
include 'plot3d_operator.f90'
include 'CUBE_mod.f90'

program CUBE2USG
    use CUBE_mod
    use vtkMesh_operator_m
    implicit none
    character(50) F_fname, USG_fname
    integer i, CUBEID, nodeID(3), n, num_node
    real X(3)

    print *, 'F_FileName ?'
    read(5,*) F_fname

    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call read_CUBE_data(F_fname, '')

    call read_VTK_mesh(USG_fname, meshONLY=.true.)

    do i = 0, size(cell_array)-1
        X(:) = 0.0
        num_node = size(cell_array(i)%nodeID)
        do n = 1, num_node
            X(:) = X(:) + node_array(cell_array(i)%nodeID(n))%coordinate(:)
        end do
        X(:) = X(:) / real(num_node)
        CUBEID = get_cube_contains(X)    
        nodeID(:) = nearest_node(X, CUBEID)
        cell_array(i)%vector(:) = get_velocity_f(nodeID, CUBEID)
    end do

    i = len_trim(F_fname)
    call output_VTK_mesh(F_fname(:i-2)//'.vtk')
    
end program CUBE2USG