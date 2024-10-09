module path_operator_m
    implicit none

    ! character(7), parameter, private :: OS = 'Linux'

    contains

    subroutine make_directory(path)
        character(*), intent(in) :: path
        character(:), allocatable :: directory
        character(len=7) OS

        call get_os(OS)
    
        select case(trim(OS))
            case ('Linux')  !for_Linux
                directory = replace_str(path, from='\', to='/')
                call system('mkdir -p "' // directory // '"')

            case ('Windows')  !for_Windows
                directory = replace_str(path, from='/', to='\')
                call system('md "' // directory // '"')

            case default
                print*, 'OS ERROR : ', OS
                error stop
                
        end select

    end subroutine make_directory

    
    subroutine get_os(OS)
        character(len=7), intent(out) :: OS
        character(len=128) env_os

        call getenv('OS', env_os)
        env_os = adjustl(trim(env_os))

        if(index(env_os, 'Windows') > 0) then
            OS = 'Windows'
        else
            OS = 'Linux'
        end if

    end subroutine get_os
    
    subroutine get_DirFromPath(path, directory, filename)
        character(*), intent(in) :: path
        character(:), intent(out), allocatable :: directory
        character(:), intent(out), allocatable , optional :: filename
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

        if(present(filename)) filename = trim(path(i+1:))
        directory = path(:i)

        print*, 'Path= ', trim(path)
        print*, 'Directory= ', directory
        if(present(filename)) print*, 'Filename= ', filename

    end subroutine

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

    end function

end module path_operator_m