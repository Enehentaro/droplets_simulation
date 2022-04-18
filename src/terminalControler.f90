module terminalControler_m
    implicit none
    private

    character(1), parameter :: esc = achar(27)
    integer, parameter :: stdOut = 6
    character(:), allocatable :: format_str

    public reset_formatTC, set_formatTC, print_sameLine

    interface print_sameLine
        module procedure print_sameLine_int, print_sameLine_real
    end interface

    !============================================================================
    ! How to use EXAMPLE:

    ! call set_formatTC('("CHECK halfFace [ #group : ",i6," / ",i6," ]")')
    ! do groupID = 1, num_group
    !     call print_sameLine([groupID, num_group])
    ! end do

    !============================================================================
            
    contains

    subroutine reset_formatTC
        
        format_str = ''

    end subroutine 

    subroutine set_formatTC(fmt_str)
        character(*), intent(in) :: fmt_str

        if(fmt_str == format_str) return

        format_str = fmt_str
        write(stdOut, '()')     !改行（あとで戻ってくるため）

    end subroutine 

    subroutine print_sameLine_int(array)
        integer, intent(in) :: array(:)

        write(stdOut, "(a)", advance='no') esc//'M'     !カーソルが1行戻る
        write(stdOut, format_str) array(:)

    end subroutine

    subroutine print_sameLine_real(array)
        real, intent(in) :: array(:)

        write(stdOut, "(a)", advance='no') esc//'M'     !カーソルが1行戻る
        write(stdOut, format_str) array(:)

    end subroutine

end module
