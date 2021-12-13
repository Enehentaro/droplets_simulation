module boxCounter_m
    implicit none

    type boxCounter
        real center(3), width(3), min_cdn(3), max_cdn(3)
        logical, allocatable :: dropletFlag(:)

        contains

        procedure add_dropletFlag
        procedure get_id_array
    end type

    contains

    function get_box_array(dir, num_drop) result(new_box_array)
        use simpleFile_reader
        type(boxCounter), allocatable :: new_box_array(:)
        character(*), intent(in) :: dir
        integer, intent(in) :: num_drop
        double precision, allocatable :: boxSize_mat(:,:)
        integer i, num_box

        call read_CSV(filename=dir//'/boxList.csv', matrix=boxSize_mat, header=.true.)

        num_box = size(boxSize_mat, dim=2)

        allocate(new_box_array(num_box))

        do i = 1, num_box
            new_box_array(i)%center(:) = real(boxSize_mat(1:3, i))
            new_box_array(i)%width(:) = real(boxSize_mat(4:6, i))
            new_box_array(i)%min_cdn(:) = new_box_array(i)%center(:) - new_box_array(i)%width(:)*0.5
            new_box_array(i)%max_cdn(:) = new_box_array(i)%center(:) + new_box_array(i)%width(:)*0.5

            if((new_box_array(i)%width(1) <= 0.d0)&
                .or.(new_box_array(i)%width(2) <= 0.d0).or.(new_box_array(i)%width(3) <= 0.d0)) then
                    print*, 'ERROR boxSize :', new_box_array(i)%center(:), new_box_array(i)%width(:)
                    stop
            end if

            allocate(new_box_array(i)%dropletFlag(num_drop), source=.false.)
        end do

    end function

    subroutine add_dropletFlag(self, id_array)
        class(boxCounter) self
        integer, intent(in) :: id_array(:)

        self%dropletFlag(id_array) = .true.

    end subroutine

    function get_id_array(self) result(id_array)
        class(boxCounter) self
        integer, allocatable :: id_array(:)
        integer i, cnt, n_size

        n_size = count(self%dropletFlag)
        allocate(id_array(n_size))

        cnt = 0
        do i = 1, size(self%dropletFlag)
            if(self%dropletFlag(i)) then
                cnt = cnt + 1
                id_array(cnt) = i
            end if
        end do

    end function

end module boxCounter_m