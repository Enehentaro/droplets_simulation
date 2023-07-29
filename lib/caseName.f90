module caseName_m
    implicit none

    contains

    subroutine case_check(caseName_array)
        !!case名をキーボードから取得する。
        !!TXTファイルを指定すると、それを全行読み込んで配列に格納。
        use simpleFile_reader
        character(*), allocatable, intent(out) :: caseName_array(:)
        character(255) caseName
        integer i
        logical existance
        character(21), parameter :: conditionFName = 'condition.nml'

        print*, 'Case Name ?'
        read(5, '(A)') caseName
        if(index(caseName, '.txt') > 0) then
            call read_textRecord(trim(caseName), caseName_array)
        else
            caseName_array = [caseName]
        end if

        do i = 1, size(caseName_array)
            inquire(file=trim(caseName_array(i))//'/'//conditionFName, exist=existance)
            if(.not.existance) then
                print*, 'Case:[ ', trim(caseName_array(i)), ' ] is not found.'
                error stop
            end if
        end do

    end subroutine

end module caseName_m
