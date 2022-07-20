module caseName_m
    implicit none
    integer, target :: nowCase
    character(30), allocatable, private :: case_array(:)

    contains

    subroutine case_check(num_case)
        use simpleFile_reader
        use filename_m, only : conditionFName
        character(30) caseName
        integer, intent(out) :: num_case
        integer i
        logical existance

        print*, 'Case Name ?'
        read(5, '(A)') caseName
        if(index(caseName, '.txt') > 0) then
            call read_textRecord(trim(caseName), case_array)
            num_case = size(case_array)
        else
            num_case = 1
            allocate(case_array(1))
            case_array(1) = caseName
        end if

        do i = 1, num_case
            inquire(file=trim(case_array(i))//'/'//conditionFName, exist=existance)
            if(.not.existance) then
                print*, 'Case:[ ', trim(case_array(i)), ' ] is not found.'
                error stop
            end if
        end do

    end subroutine

    function get_caseName() result(caseName)
        character(:), allocatable :: caseName

        caseName = trim(case_array(nowCase))

        print*, '#', nowCase, '[',caseName,']'

    end function

end module caseName_m
