module cases_reader
    implicit none
    
    type :: case_condition
        character(99) :: path, path2
        integer T, RH
    end type case_condition
    type(case_condition), allocatable, private :: case_matrix(:)

    integer, private :: num_cases = 1

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
            
        end if

    end function check_cases

    subroutine read_cases(FNAME)
        character(*), intent(in) :: FNAME
        ! character(99) A
        integer n_unit, i, ios

        print*, 'READ:', FNAME
            
        open (newunit=n_unit, file=FNAME, status='old')
    
            num_cases = 0
            read(n_unit, '()') !ヘッダーの読み飛ばし
            do        !レコード数を調べるループ
                read(n_unit, '()', iostat=ios)
                if(ios/=0) exit
                num_cases = num_cases + 1
            end do

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

        print*, path, directory, filename

    end subroutine set_dir_from_path

    subroutine replace_str( str, from, to )
        character (*),intent (inout) :: str
        character (1),intent (in) :: from, to
        integer :: i, l

        l = len_trim(str)
        do i=1, l
            if ( str(i:i) == from ) str(i:i) = to
        end do

    end subroutine replace_str

    integer function get_temperature(index)
        integer, intent(in) :: index

        get_temperature = case_matrix(index)%T

    end function get_temperature

    integer function get_humidity(index)
        integer, intent(in) :: index

        get_humidity = case_matrix(index)%RH

    end function get_humidity

    subroutine set_case_path(case_path, index)
        character(*), intent(inout) :: case_path
        integer, intent(in) :: index

        case_path = trim(case_matrix(index)%path)

    end subroutine set_case_path

    subroutine set_case_path2(case_path, index)
        character(*), intent(inout) :: case_path
        integer, intent(in) :: index

        case_path = trim(case_matrix(index)%path2)

    end subroutine set_case_path2

end module cases_reader