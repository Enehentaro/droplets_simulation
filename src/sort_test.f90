program sort_test
    use sort_m
    implicit none
    integer, parameter :: arraySize = 100
    type(content_t) array(arraySize), array_sorted(arraySize)
    integer i
    real ct1, ct2

    call random_number(array(:)%value)
    do i = 1, arraySize
        array(i)%originID = i
    end do

    call cpu_time(ct1)
    call heap_sort(array, array_sorted)
    call cpu_time(ct2)

    do i = 1, arraySize
        print '(2(i6, x, f20.16, 4x, "|"))', array(i)%originID, array(i)%value, array_sorted(i)%originID, array_sorted(i)%value
    end do

    print*, ct2 - ct1, 'sec'
    
end program sort_test