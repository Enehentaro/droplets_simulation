module cases_reader
    implicit none
    
    type :: case_condition
        character(99) :: path, path2
        integer T, RH
    end type case_condition
    type(case_condition), allocatable, private :: case_matrix(:)

    integer, private :: num_cases = 1

    logical :: cases_read_flag = .false.

    contains

    integer function check_cases(fname)
        character(*), intent(in) :: fname
        integer i

        i = index(fname, '.csv')

        if(i <= 0) then
            print*, 'Normal_program'
            check_cases = 1

        else
            call read_cases(fname)
            check_cases = num_cases
            cases_read_flag = .true.
            
        end if

    end function check_cases

    subroutine read_cases(FNAME)
        use csv_reader
        character(*), intent(in) :: FNAME
        integer n_unit, i
        integer :: mat_size(2)

        print*, 'READ:', FNAME
            
        open (newunit=n_unit, file=FNAME, status='old')
            
            mat_size = get_size(n_unit, header_flag=.true.)

            num_cases = mat_size(2)

            allocate(case_matrix(num_cases))

            rewind(n_unit)  ! ファイルの最初に戻る
            print *, 'NumCases =', num_cases
            read(n_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, num_cases        !本読み込み
                read(n_unit, *) case_matrix(i)%path, case_matrix(i)%path2, case_matrix(i)%T, case_matrix(i)%RH
                print *, trim(case_matrix(i)%path), ' ', trim(case_matrix(i)%path2), case_matrix(i)%T, case_matrix(i)%RH
            end do

        close (n_unit)


    end subroutine read_cases

    subroutine set_dir_from_path(path, directory, filename)
        character(*), intent(in) :: path
        character(*), intent(inout) :: directory
        character(*), intent(inout), optional :: filename
        character(1) delimiter
        integer i

        if(index(path, '/') > 0) then
            delimiter = '/'

        else if(index(path, '\') > 0) then
            delimiter = '\'

        else
            print*, 'Delimiter was not found.'
            if(present(filename)) filename = path
            directory = ''
            return

        end if

        i = index(path, delimiter, back=.true.)

        if(present(filename)) filename = path(i+1:)
        directory = path(:i)

        print*, 'Path= ', trim(path)
        print*, 'Directory= ', trim(directory)
        print*, 'Filename= ', trim(filename)

    end subroutine set_dir_from_path

    function replace_str( str, from, to )
        character (*),intent (in) :: str
        character (1),intent (in) :: from, to
        character (len_trim(str)) :: replace_str
        integer :: i, l

        replace_str = str
        l = len_trim(str)
        do i = 1, l
            if ( str(i:i) == from ) replace_str(i:i) = to
        end do

    end function replace_str

    integer function get_temperature(index)
        integer, intent(in) :: index

        get_temperature = case_matrix(index)%T

    end function get_temperature

    integer function get_humidity(index)
        integer, intent(in) :: index

        get_humidity = case_matrix(index)%RH

    end function get_humidity

    function get_case_path(index)
        integer, intent(in) :: index
        character(len_trim(case_matrix(index)%path)) get_case_path

        get_case_path = trim(case_matrix(index)%path)

    end function get_case_path

    function get_case_path2(index)
        integer, intent(in) :: index
        character(len_trim(case_matrix(index)%path2)) get_case_path2

        get_case_path2 = trim(case_matrix(index)%path2)

    end function get_case_path2

end module cases_reader