program main
    use sort_m

    implicit none
    integer n_unit
    integer i
    character(50) txtFName, outputFName
    real, allocatable :: sample(:),after_sample(:)

    txtFName = "/home/master/Documents/kd-tree-git/droplets_simulation/Test/sample.txt"
    outputFName = "/home/master/Documents/kd-tree-git/droplets_simulation/Test/output.txt"
    allocate(sample(10))

    open(newunit = n_unit, file = txtFName, status = 'old')
        do i = 1, 10
            read(n_unit,*) sample(i)
        end do
    close(n_unit)

    call construct_heap(sample, after_sample)

    open(newunit = n_unit, file = outputFName, status = 'replace')
        do i = 1, 10
            write(n_unit,*) after_sample(i)
        end do
    close(n_unit)

end program main