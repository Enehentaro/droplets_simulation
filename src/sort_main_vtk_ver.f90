program sortMain_vtk_ver

    use sort_m
    use kdTree_m
    use unstructuredGrid_mod
    implicit none 

    type(UnstructuredGrid) grid
    type(content_t), allocatable :: before(:), after(:)
    type(content_t), allocatable :: leftChild(:), rightChild(:), kdTree_cell(:)

    integer :: n_unit 
    integer :: i, iimx, k, kkmx
    character(50) vtkFName

    vtkFName = "Test/sample.vtk"
            
    grid = UnstructuredGrid_(vtkFName)

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(before(size(grid%CElls)))
    allocate(after(size(grid%CElls)))

    do i = 1, iimx
        before(i)%originID = i
        before(i)%value = grid%CELLs(i)%center(1)
    end do

    open(newunit = n_unit, file = "Test/No1_before.txt", status = 'replace')
        do i = 1, iimx
            write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(f12.5)') before(i)%value
        end do
    close(n_unit)

    call heap_sort(before, after)

    open(newunit = n_unit, file = "Test/No2_after.txt", status = 'replace')
        do i = 1, iimx 
            write(n_unit,'(I0)', advance='no') after(i)%originID
            write(n_unit,'(f12.5)') after(i)%value
        end do
    close(n_unit)

    call cellDivider(grid, after, leftChild, kdTree_cell, rightChild)

    open(newunit = n_unit, file = "Test/No3_divide_after.txt", status = 'replace')
        do i = 1, size(leftChild) 
            write(n_unit,'(I0)', advance='no') leftChild(i)%originID
            write(n_unit,'(f12.5)') leftChild(i)%value
        end do
        write(n_unit,'(A)')
        write(n_unit,'(I0)', advance='no') kdTree_cell(1)%originID
        write(n_unit,'(f12.5)') kdTree_cell(1)%value
        write(n_unit,'(A)')
        do i = 1, size(rightChild) 
            write(n_unit,'(I0)', advance='no') rightChild(i)%originID
            write(n_unit,'(f12.5)') rightChild(i)%value
        end do
    close(n_unit)

end program sortMain_vtk_ver