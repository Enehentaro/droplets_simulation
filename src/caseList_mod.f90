module case_list_m
    implicit none
    character(30), allocatable, private :: case_list(:)
    character(13), parameter, private :: caseList_Fname = 'case_list.txt'

    contains

    subroutine case_check(num_case)
        use filename_mod
        character(30) case_name
        integer, intent(out) :: num_case
        integer i
        logical existance

        print*, 'Case Name ?'
        read(5, *) case_name
        if(trim(case_name)==caseList_Fname) then
            call read_case_list
            num_case = size(case_list)
        else
            num_case = 1
            allocate(case_list(1))
            case_list(1) = case_name
        end if

        do i = 1, num_case
            inquire(file=trim(case_list(i))//'/'//conditionFName, exist=existance)
            if(.not.existance) then
                print*, 'Case : ', trim(case_list(i)), ' is not found.'
                stop
            end if
        end do

    end subroutine case_check

    subroutine read_case_list
        integer n_unit, num_rec, ios, i
        character(10) A

        open(newunit=n_unit, file=caseList_Fname, status='old')
            num_rec = 0
            do        !レコード数を調べるループ
                read(n_unit, '(A)', iostat=ios) A !ファイル終端であればiosに-1が返る
                if((trim(A) == '').or.(ios/=0)) exit    !終端もしくは空白行であればループ脱出
                num_rec = num_rec + 1

            end do

            allocate(case_list(num_rec))

            rewind(n_unit)

            do i = 1, num_rec
                read(n_unit, '(A)', iostat=ios) case_list(i) !ファイル終端であればiosに-1が返る
            end do

        close(n_unit)
    end subroutine

    function get_case_name(nc) result(case_name)
        integer, intent(in) :: nc
        character(:), allocatable :: case_name

        case_name = trim(case_list(nc))

    end function get_case_name

end module case_list_m