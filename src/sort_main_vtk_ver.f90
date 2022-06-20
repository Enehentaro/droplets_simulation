program sortMain_vtk_ver

    use sort_m
    use kdTree_m
    use unstructuredGrid_mod
    implicit none 

    type(UnstructuredGrid) grid
    type(content_t), allocatable :: before(:), after(:)
    type(nodeTree_t), allocatable :: node_tree(:)
    real droplet_position(3)
    integer nearest_ID

    integer n_unit 
    integer i, iimx, k, kkmx
    character(50) vtkFName

    vtkFName = "Test/sample.vtk"
            
    grid = UnstructuredGrid_(vtkFName)

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(before(size(grid%CElls)))

    call system('mkdir -p -v Test_check')

    do i = 1, iimx
        before(i)%originID = i
        before(i)%coordinate(1) = grid%CELLs(i)%center(1)
        before(i)%coordinate(2) = grid%CELLs(i)%center(2)
        before(i)%coordinate(3) = grid%CELLs(i)%center(3)
    end do

    open(newunit = n_unit, file = "Test_check/before.txt", status = 'replace')
        do i = 1, iimx
            write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(3(f12.5))') before(i)%coordinate(1), before(i)%coordinate(2), before(i)%coordinate(3)
        end do
    close(n_unit)

    call create_kdtree(before, node_tree)

    do i = 1, 113
        print*, i
        droplet_position(:) = before(node_tree(i)%cell_ID)%coordinate(:)
        call search_kdtree(before, node_tree, droplet_position, nearest_ID)
    end do

end program sortMain_vtk_ver