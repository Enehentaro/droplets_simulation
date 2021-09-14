module cases_reader
    implicit none
    
    type :: case_condition
        character :: path*99, fname*20
        integer T, RH
    end type case_condition
    type(case_condition), allocatable, private :: case_matrix(:)

    integer, private :: num_cases = 1

    contains

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
                ! read(n_unit, '(A)') A
                ! call replace_str( A, '/', '\' )
                ! read(A, *) case_matrix(i)%path, case_matrix(i)%fname, case_matrix(i)%T, case_matrix(i)%RH
                read(n_unit, *) case_matrix(i)%path, case_matrix(i)%fname, case_matrix(i)%T, case_matrix(i)%RH
                ! if(trim(OS)=='Linux') call replace_str( case_matrix(i)%path, '\', '/' )
                print *, trim(case_matrix(i)%path), case_matrix(i)%fname, case_matrix(i)%T, case_matrix(i)%RH
            end do

        close (n_unit)


    end subroutine read_cases

    subroutine replace_str( str, from, to )
        character (*),intent (inout) :: str
        character (1),intent (in) :: from, to
        integer :: i, l

        l = len_trim(str)
        do i=1, l
              if ( str(i:i) == from ) str(i:i) = to
        end do

    end subroutine replace_str

    integer function get_num_cases()

        get_num_cases = num_cases

    end function get_num_cases

    integer function get_temperature(index)
        integer, intent(in) :: index

        get_temperature = case_matrix(index)%T

    end function get_temperature

    integer function get_humidity(index)
        integer, intent(in) :: index

        get_humidity = case_matrix(index)%RH

    end function get_humidity

    subroutine set_path_out_base(path_out_base, index)
        character(*), intent(inout) :: path_out_base
        integer, intent(in) :: index
        integer i

        path_out_base = trim(case_matrix(index)%fname)
        i = len_trim(path_out_base)
        if(path_out_base(i:i) == '\') path_out_base(i:i) = ' '
        path_out_base = '..\' // trim(path_out_base) // '_virus\'

    end subroutine set_path_out_base

    subroutine set_head_out(head_out, index)
        character(*), intent(inout) :: head_out
        integer, intent(in) :: index

        head_out = trim(case_matrix(index)%fname)

    end subroutine set_head_out

    subroutine set_PATH_AIR(PATH_AIR, index)
        character(*), intent(inout) :: PATH_AIR
        integer, intent(in) :: index

        PATH_AIR = trim(case_matrix(index)%path)

    end subroutine set_PATH_AIR

    subroutine set_HEAD_AIR(HEAD_AIR, index)
        character(*), intent(inout) :: HEAD_AIR
        integer, intent(in) :: index

        HEAD_AIR = trim(case_matrix(index)%fname)

    end subroutine set_HEAD_AIR

end module cases_reader