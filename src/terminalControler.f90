module terminalControler_m
    implicit none
    character(1), parameter, private :: esc = achar(27)
    integer, parameter, private :: stdOut = 6
    character(:), allocatable, private :: format_str
            
    contains

    subroutine set_formatTC(fmt_str)
        character(*), intent(in) :: fmt_str

        if(fmt_str == format_str) return

        format_str = fmt_str
        write(stdOut, '()')

    end subroutine 

    subroutine print_sameLine(int_array)
        integer, intent(in) :: int_array(:)

        write(stdOut, "(a)", advance='no') esc//'M'
        write(stdOut, format_str) int_array(:)

    end subroutine
end module