module simpleFile_reader
    implicit none
    private
    !簡単なファイル（CSV、TXTなど）の操作手続き集モジュール

    interface read_CSV
        module procedure read_csv_dble, read_csv_int, read_csv_char
    end interface

    public read_CSV, read_textRecord

    contains

    subroutine read_csv_char(filename, matrix, column, header)
        integer i, Num_unit
        character(*), intent(in) :: filename
        character(*), intent(inout), allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old', action='read')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                ! print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine read_csv_char

    subroutine read_csv_dble(filename, matrix, column, header)
        integer i, Num_unit
        character(*), intent(in) :: filename
        double precision, intent(inout), allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old', action='read')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                ! print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine read_csv_dble

    subroutine read_csv_int(filename, matrix, column, header)
        integer i, Num_unit
        character(*), intent(in) :: filename
        integer, intent(inout), allocatable :: matrix(:,:)
        integer :: mat_size(2)
        integer, intent(in), optional :: column
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old', action='read')

            mat_size = get_size_csv(Num_unit, header_flag)

            if(present(column)) mat_size(1) = column

            allocate(matrix(mat_size(1), mat_size(2)))
            print *, 'Size =', mat_size(:)

            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, mat_size(2)        !本読み込み
                read (Num_unit, *) matrix(:,i)
                ! print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine read_csv_int

    subroutine read_textRecord(filename, array)
        integer i, Num_unit
        character(*), intent(in) :: filename
        character(*), intent(out), allocatable :: array(:)
        integer :: num_record

        print*, 'simpleREADER : ', filename

        open (newunit=Num_unit, file=filename, status='old', action='read')

            num_record = get_num_records(Num_unit, header_flag=.false.)

            allocate(array(num_record))
            print *, '#Records =', num_record

            do i = 1, num_record
                read (Num_unit, '(A)') array(i)
            end do

        close (Num_unit)

    end subroutine

    function get_size_csv(Num_unit, header_flag) result(mat_size)
        integer, intent(in) :: Num_unit
        logical, intent(in) :: header_flag
        integer :: mat_size(2)
        character(255) A

        mat_size(:) = 0

        read(Num_unit, '(A)') A
        mat_size(1) = get_num_columns(A)

        mat_size(2) = get_num_records(Num_unit, header_flag)

    end function get_size_csv

    integer function get_num_records(Num_unit, header_flag)
        integer, intent(in) :: Num_unit
        logical, intent(in) :: header_flag
        character(20) A
        integer ios

        rewind(Num_unit)  ! ファイルの最初に戻る
        get_num_records = 0
        if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
        do        !レコード数を調べるループ
            read(Num_unit, '(A)', iostat=ios) A !ファイル終端であればiosに-1が返る
            if((trim(A) == '').or.(ios/=0)) exit    !終端もしくは空白行であればループ脱出
            get_num_records = get_num_records + 1
        end do

        rewind(Num_unit)  ! ファイルの最初に戻る

    end function

    integer function get_num_columns(str)
        character(*), intent(in) :: str
        character(1) delimiter
        integer i, j

        delimiter = ','
        if(index(str, ',') == 0) delimiter = ' '

        get_num_columns = 1
        i = 1
        do
            j = index(str(i:), delimiter)
            if(j > 0) then
                if(j > 1) get_num_columns = get_num_columns + 1
                i = i + j
            else
                exit
            end if

        end do
        
    end function

end module simpleFile_reader