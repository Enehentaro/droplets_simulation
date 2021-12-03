module caseNameList_m
    implicit none
    integer, target :: nowCase
    character(30), allocatable, private :: case_list(:)
    ! character(13), parameter, private :: caseList_Fname = 'case_list.txt'

    contains

    subroutine case_check(num_case)
        use filename_mod
        character(30) caseName
        integer, intent(out) :: num_case
        integer i
        logical existance

        print*, 'Case Name ?'
        read(5, '(A)') caseName
        inquire(file=trim(caseName), exist=existance)
        if(existance) then
            call read_case_list(trim(caseName))
            num_case = size(case_list)
        else
            num_case = 1
            allocate(case_list(1))
            case_list(1) = caseName
        end if

        do i = 1, num_case
            inquire(file=trim(case_list(i))//'/'//conditionFName, exist=existance)
            if(.not.existance) then
                print*, 'Case:[ ', trim(case_list(i)), ' ] is not found.'
                stop
            end if
        end do

    end subroutine

    subroutine read_case_list(fname)
        integer n_unit, num_rec, ios, i
        character(*), intent(in) :: fname
        character(10) A

        print*, 'READ: ', fname
        open(newunit=n_unit, file=fname, status='old')
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

    function get_caseName() result(caseName)
        character(:), allocatable :: caseName

        caseName = trim(case_list(nowCase))

    end function

end module caseNameList_m