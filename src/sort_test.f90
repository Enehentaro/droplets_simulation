program sort_test
    use sort_m
    implicit none
    integer, parameter :: arraySize = 10
    type(content_t) array(arraySize), array_sorted(arraySize)
    integer i

    call random_number(array(:)%value)
    do i = 1, arraySize
        array(i)%originID = i
    end do

    do i = 1, arraySize
        print*, array(i)%originID, array(i)%value
    end do

    print '("==================================================")'

    call heap_sort(array, array_sorted)

    do i = 1, arraySize
        print*, array_sorted(i)%originID, array_sorted(i)%value
    end do
    
end program sort_test