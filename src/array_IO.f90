module array_IO_m

    contains

    subroutine output_array_asBinary(fname, array)
        character(*), intent(in)  :: fname
        real, intent(in) :: array(:,:)
        integer n_unit

        print*, 'output_array : ', fname

        open(newunit=n_unit, file=fname, form='unformatted', status='replace')

            write(n_unit) shape(array)

            write(n_unit) array

        close(n_unit)

    end subroutine

    subroutine read_array_asBinary(fname, array)
        character(*), intent(in)  :: fname
        real, allocatable, intent(out) :: array(:,:)
        integer n_unit, arrayShape(2)

        print*, 'read_array : ', fname

        open(newunit=n_unit, file=fname, form='unformatted', status='old')

            read(n_unit) arrayShape(:)

            allocate(array(arrayShape(1), arrayShape(2)))

            read(n_unit) array

        close(n_unit)

    end subroutine

end module