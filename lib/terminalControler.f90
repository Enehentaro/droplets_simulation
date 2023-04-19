module terminalControler_m
    !!ターミナル操作を扱うモジュール
    implicit none
    private

    character(1), parameter :: esc = achar(27)
    integer, parameter :: stdOut = 6
    character(:), allocatable :: format_str

    public reset_formatTC, set_formatTC, print_progress, print_sameLine

    interface print_progress
        !!進捗を表示する
        module procedure print_progress_int, print_progress_real
    end interface

    !|##Example
    !```Fortran
    !program sample
    !   use terminalControler_m
    !
    !   call set_formatTC('("PROGRESS [ #progress : ",i6," / ",i6," ]")')
    !   do i = 1, imax
    !       call print_progress([i, imax])
    !   end do
    !
    !end program sample
    !```

    contains

    subroutine reset_formatTC
        
        format_str = ''

    end subroutine 

    subroutine set_formatTC(fmt_str)
        !!進捗を表示するためのフォーマットを指定
        !!指定時に改行が起こる（あとで戻ってくるため）
        character(*), intent(in) :: fmt_str

        if(allocated(format_str)) then
            if(fmt_str == format_str) return
        end if

        format_str = fmt_str
        write(stdOut, '()')     !改行（あとで戻ってくるため）

    end subroutine 

    subroutine print_progress_int(array)
        integer, intent(in) :: array(:)

        write(stdOut, "(a)", advance='no') esc//'M'     !カーソルが1行戻る
        write(stdOut, format_str) array(:)

    end subroutine

    subroutine print_progress_real(array)
        real, intent(in) :: array(:)

        write(stdOut, "(a)", advance='no') esc//'M'     !カーソルが1行戻る
        write(stdOut, format_str) array(:)

    end subroutine

    subroutine print_sameLine(str)
        character(*), intent(in) :: str

        write(stdOut, "(a)", advance='no') esc//'M'     !カーソルが1行戻る
        write(stdOut, "(a)") str

    end subroutine

end module
