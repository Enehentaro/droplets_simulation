program sortMain_vtk_ver
    use kdTree_m
    use unstructuredGrid_mod
    implicit none 

    type(UnstructuredGrid) grid
    real, allocatable :: xyz(:,:)
    type(node_in_kdTree_t), allocatable :: kdTree(:)
    real droplet_position(3)
    integer nearest_ID

    integer n_unit 
    integer i, iimx, kkmx
    character(50) vtkFName

    vtkFName = "Test/sample2.vtk"
            
    grid = UnstructuredGrid_(vtkFName)

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(xyz(3, size(grid%CElls)))

    call system('mkdir -p -v Test_check')

    do i = 1, iimx
        xyz(:, i) = grid%CELLs(i)%center(:)
    end do

    open(newunit = n_unit, file = "Test_check/before.txt", status = 'replace')
        do i = 1, iimx
            ! write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(3(f12.5))') xyz(:, i)
        end do
    close(n_unit)

    call create_kdtree(xyz, kdTree)

    do i = 1, 63
        print*, i
        droplet_position(:) = xyz(:, i)
        call search_kdtree(xyz, kdTree, droplet_position, nearest_ID)
    end do

end program sortMain_vtk_ver