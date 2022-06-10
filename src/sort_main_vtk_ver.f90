program sortMain_vtk_ver

    use sort_m
    use unstructuredGrid_mod
    implicit none 

    type(UnstructuredGrid) grid
    
    type(content), allocatable :: after_sample(:)

    integer :: n_unit 
    integer :: i, iimx, k, kkmx
    character(50) vtkFName, before_outputFName, after_outputFName
    real, allocatable :: sample(:)
    integer, allocatable :: preCellID(:)

    vtkFName = "Test/sample.vtk"
    before_outputFName = "Test/before_output.txt"
    after_outputFName = "Test/after_output.txt"
            
    call grid%read_VTK(vtkFName,meshOnly=.true.)
    call grid%set_gravity_center()

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(preCellID(iimx))

    open(newunit = n_unit, file = before_outputFName, status = 'replace')
        do i = 1, iimx
            preCellID(i) = i
            write(n_unit,'(I3)', advance='no') preCellID(i)
            write(n_unit,'(f12.5)') grid%CELLs(i)%center(1)
        end do
    close(n_unit)

    allocate(sample(size(grid%CElls)))
    sample(:) = grid%CElls(:)%center(1) 

    call heap_sort(sample, preCellID, after_sample)

    open(newunit = n_unit, file = after_outputFName, status = 'replace')
        do i = 1, iimx 
            write(n_unit,'(I0)', advance='no') after_sample(i)%originID
            write(n_unit,'(f12.5)') after_sample(i)%axis
        end do
    close(n_unit)

end program sortMain_vtk_ver