module boxCounter_m
    implicit none

    type box_t
        real center(3), width(3), min_cdn(3), max_cdn(3)
    end type

    type, extends(box_t) :: boxCounter
        logical, allocatable :: Flag(:)
        contains
        procedure add_Flag
        procedure get_FlagID
    end type

    contains

    function get_boxArray_fromCSV(fname) result(new_box_array)
        use simpleFile_reader
        type(box_t), allocatable :: new_box_array(:)
        character(*), intent(in) :: fname
        double precision, allocatable :: boxSize_mat(:,:)
        integer i, num_box

        call read_CSV(filename=fname, matrix=boxSize_mat, header=.true.)

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

        end do

    end function

    function get_boxCounterArray_fromCSV(fname, num_Flag) result(boxCounter_array)
        type(box_t), allocatable :: box_array(:)
        type(boxCounter), allocatable :: boxCounter_array(:)
        character(*), intent(in) :: fname
        integer, intent(in) :: num_Flag
        integer i, num_box

        box_array =  get_boxArray_fromCSV(fname)
        num_box = size(box_array)

        do i = 1, num_box
            boxCounter_array(i)%box_t = box_array(i)
            allocate(boxCounter_array(i)%Flag(num_Flag), source=.false.)
        end do

    end function

    subroutine add_Flag(self, id_array)
        class(boxCounter) self
        integer, intent(in) :: id_array(:)

        self%Flag(id_array) = .true.

    end subroutine

    function get_FlagID(self) result(id_array)
        class(boxCounter) self
        integer, allocatable :: id_array(:)
        integer i, cnt, n_size

        n_size = count(self%Flag)
        allocate(id_array(n_size))

        cnt = 0
        do i = 1, size(self%Flag)
            if(self%Flag(i)) then
                cnt = cnt + 1
                id_array(cnt) = i
            end if
        end do

    end function

end module boxCounter_m
