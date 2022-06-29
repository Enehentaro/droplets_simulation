program sort_test
    use sort_m
    use array_m
    implicit none
    real, parameter :: test1(10) = [-9.0, -1.0, 0.0, 3.0, 5.0, 7.2, 14.4, 99.9, 122.5, 255.0]
    real, parameter :: test2(11) = [-99.0, -9.0, -1.0, 0.0, 3.0, 5.0, 7.2, 14.4, 99.9, 122.5, 255.0]

    call testing(test1)
    print '("==============================================================")'
    call testing(test2)

    contains

    subroutine testing(array_correct)
        real, intent(in) :: array_correct(:)
        real tmp(size(array_correct))
        type(content_t) array(size(array_correct)), array_sorted(size(array_correct))
        integer i

        tmp = FisherYates_shuffle(array_correct)
        ! print *, tmp
        array = real2content(tmp)
        call heap_sort(array, array_sorted)

        ! if(.not.isEqual(array_sorted(:)%value, array_correct)) error stop
        if(.not.all(array_sorted(:)%value == array_correct)) error stop

        do i = 1, size(array_sorted)
            print '(2(i6, x, f20.16, 4x, "|"))', array(i)%originID, array(i)%value, array_sorted(i)%originID, array_sorted(i)%value
        end do

    end subroutine
    
end program sort_test