module csv_reader
    implicit none

    contains

    function read_csv_dble(filename, column, header) result(matrix)
        integer i, Num_unit
        character(*), intent(in) :: filename
        double precision, allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            rewind (Num_unit)  ! ファイルの最初に戻る

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end function read_csv_dble

    subroutine csv_reader_dble(filename, matrix, column, header)
        integer i, Num_unit
        character(*), intent(in) :: filename
        double precision, intent(inout), allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            rewind (Num_unit)  ! ファイルの最初に戻る

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine csv_reader_dble

    subroutine csv_reader_int(filename, matrix, column, header)
        integer i, Num_unit
        character(*), intent(in) :: filename
        integer, intent(inout), allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            rewind (Num_unit)  ! ファイルの最初に戻る

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine csv_reader_int

    function get_size_csv(Num_unit, header_flag) result(mat_size)
        integer, intent(in) :: Num_unit
        logical, intent(in) :: header_flag
        integer :: mat_size(2)
        character(99) A
        integer ios

        mat_size(:) = 0
        if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
        do        !レコード数を調べるループ
            read(Num_unit, '(A)', iostat=ios) A !ファイル終端であればiosに-1が返る
            if((trim(A) == '').or.(ios/=0)) exit    !終端もしくは空白行であればループ脱出
            mat_size(2) = mat_size(2) + 1
            if(mat_size(1)==0) mat_size(1) = get_columns(A)
        end do

    end function get_size_csv

    integer function get_columns(str)
        character(*), intent(in) :: str
        character(1) :: delimiter = ','
        integer i, j

        ! if(index(str, ' ') > 0) delimiter = ' '

        get_columns = 1
        i = 1
        do
            j = index(str(i:), delimiter)
            if(j > 0) then
                if(j > 1) get_columns = get_columns + 1
                i = i + j
            else
                exit
            end if

        end do
        
    end function get_columns

end module csv_reader