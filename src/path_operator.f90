module path_operator_m
    implicit none

    ! type path_set_t
    !     character(99) case_name
    !     character(99) path2FlowDIR
    !     character(30) FlowFileName
    ! end type path_set_t

    ! type(path_set_t), allocatable :: path_list(:)

    contains

    ! subroutine check_path_list(num_case)
    !     use csv_reader
    !     use filename_mod
    !     character(99), allocatable :: path_mat(:,:)
    !     character(:), allocatable :: case_name
    !     character DIR*99, FNAME*30
    !     integer, intent(out) :: num_case
    !     integer i, n_unit

    !     print*, 'Case Name ?'
    !     read(5,*) case_name

    !     ! call read_CSV('path_list.csv', matrix=path_mat)
    !     ! num_case = size(path_mat, dim=2)
    !     num_case = 1

    !     allocate(path_list(num_case))

    !     ! do i = 1, num_case
    !     !     open(newunit=n_unit, file=path_mat(1,i), status='old')  !確認用open
    !     !     close(n_unit)
    !     !     call set_dir_from_path(path_mat(1,i), DIR, FNAME)
    !     !     open(newunit=n_unit, file=trim(DIR)//conditionFName, status='old')  !確認用open
    !     !     close(n_unit)
    !     !     open(newunit=n_unit, file=trim(DIR)//IniPositionFName, status='old')  !確認用open
    !     !     close(n_unit)
    !     !     path_list(i)%path2FlowDIR = DIR
    !     !     path_list(i)%FlowFileName = FNAME
    !     !     path_list(i)%path2outputDIR = path_mat(2,i)
    !     ! end do

    !     print*,'The Number of Cases =', num_case

    ! end subroutine check_path_list

    
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

        print*, 'Path= ', trim(path)
        print*, 'Directory= ', trim(directory)
        print*, 'Filename= ', trim(filename)

    end subroutine set_dir_from_path

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

    end function replace_str

end module path_operator_m