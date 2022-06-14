program sortMain_vtk_ver

    use sort_m
    use kdTree_m
    use unstructuredGrid_mod
    implicit none 

    type(UnstructuredGrid) grid
    type(content_t), allocatable :: before(:), after(:)

    integer n_unit 
    integer i, iimx, k, kkmx
    character(50) vtkFName

    vtkFName = "Test/sample.vtk"
            
    grid = UnstructuredGrid_(vtkFName)

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(before(size(grid%CElls)))
    allocate(after(size(grid%CElls)))

    do i = 1, iimx
        before(i)%originID = i
        before(i)%coordinate(1) = grid%CELLs(i)%center(1)
        before(i)%coordinate(2) = grid%CELLs(i)%center(2)
        before(i)%coordinate(3) = grid%CELLs(i)%center(3)
    end do

    open(newunit = n_unit, file = "Test/before.txt", status = 'replace')
        do i = 1, iimx
            write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(3(f12.5))') before(i)%coordinate(1), before(i)%coordinate(2), before(i)%coordinate(3)
        end do
    close(n_unit)

    call create_kdtree(before)

end program sortMain_vtk_ver