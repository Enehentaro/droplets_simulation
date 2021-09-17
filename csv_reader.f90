module csv_reader
    implicit none

    contains

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

            mat_size = get_size(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))

            rewind (Num_unit)  ! ファイルの最初に戻る
            print *, 'Size =', mat_size(:)
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

            mat_size = get_size(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))

            rewind (Num_unit)  ! ファイルの最初に戻る
            print *, 'Size =', mat_size(:)
            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine csv_reader_int

    function get_size(Num_unit, header_flag) result(mat_size)
        integer, intent(in) :: Num_unit
        logical, intent(in) :: header_flag
        integer :: mat_size(2)
        character(99) A
        integer ios

        mat_size(:) = 0
        if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
        do        !レコード数を調べるループ
            read(Num_unit, '(A)', iostat=ios) A !ファイル終端であればiosに-1が返る
            if((ios/=0).or.(trim(A) == '')) exit    !終端もしくは空白行があればループ脱出
            mat_size(2) = mat_size(2) + 1
            if(mat_size(1)==0) mat_size(1) = get_columns(A)
        end do

    end function get_size

    integer function get_columns(str)
        character(*), intent(in) :: str
        integer i, j

        get_columns = 1
        i = 1
        do
            j = index(str(i:), ',')
            if(j > 0) then
                get_columns = get_columns + 1
                i = i + j
            else
                exit
            end if

        end do
        
    end function get_columns

end module csv_reader