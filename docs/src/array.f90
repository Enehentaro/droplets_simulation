module array_m
    implicit none

    contains

    subroutine output_2dArray_asBinary(fname, array)
        character(*), intent(in)  :: fname
        real, intent(in) :: array(:,:)
        integer n_unit

        print*, 'output_bin2dArray : ', fname

        open(newunit=n_unit, file=fname, form='unformatted', status='replace')

            write(n_unit) shape(array)

            write(n_unit) array

        close(n_unit)

    end subroutine

    subroutine read_2dArray_asBinary(fname, array)
        character(*), intent(in)  :: fname
        real, allocatable, intent(out) :: array(:,:)
        integer n_unit, arrayShape(2)

        print*, 'read_bin2dArray : ', fname

        open(newunit=n_unit, file=fname, form='unformatted', status='old', action='read')

            read(n_unit) arrayShape(:)

            allocate(array(arrayShape(1), arrayShape(2)))

            read(n_unit) array

        close(n_unit)

    end subroutine

    subroutine read_1dArray_real(fname, array)
        character(*), intent(in)  :: fname
        real, allocatable, intent(out) :: array(:)
        integer n_unit, size

        print*, 'read_1darray : ', fname

        open(newunit=n_unit, file=fname, status='old', action='read')

            read(n_unit, *) size

            allocate(array(size))

            read(n_unit, *) array

        close(n_unit)

    end subroutine

    function mean_2dArray(array) result(mean)
        real, intent(in) :: array(:,:)
        real mean(size(array, dim=2))
        integer i

        do i = 1, size(array, dim=1)
            mean(i) = sum(array(i, :)) / size(array, dim=2)
        end do

    end function

    function FisherYates_shuffle(a) result(b)
        real, intent(in) :: a(:)
        real b(size(a)), rand, tmp
        integer i, index

        b = a

        do i = size(b), 2, -1
            call random_number(rand)
            index = int(rand * (i-1)) + 1
            ! print *, index, i

            !SWAP
            tmp = b(index)
            b(index) = b(i)
            b(i) = tmp

        end do

    end function

end module array_m