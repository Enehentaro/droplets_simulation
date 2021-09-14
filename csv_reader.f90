module csv_reader
    implicit none

    contains

    subroutine csv_reader_dble(filename, matrix, num_column, header)
        integer i, Num_unit, num_records
        integer, intent(in) :: num_column
        character(*), intent(in) :: filename
        double precision, intent(inout), allocatable :: matrix(:,:)
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old')

            num_records = get_records(Num_unit, header_flag)

            allocate(matrix(num_column, num_records))

            rewind (Num_unit)  ! ファイルの最初に戻る
            print *, 'NumRec =', num_records
            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, num_records        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine csv_reader_dble

    subroutine csv_reader_int(filename, matrix, num_column, header)
        integer i, Num_unit, num_records
        integer, intent(in) :: num_column
        character(*), intent(in) :: filename
        integer, intent(inout), allocatable :: matrix(:,:)
        logical, optional :: header
        logical :: header_flag = .true.

        if(present(header)) header_flag = header

        print*, 'CSV_READER:', filename

        open (newunit=Num_unit, file=filename, status='old')

            num_records = get_records(Num_unit, header_flag)

            allocate(matrix(num_column, num_records))

            rewind (Num_unit)  ! ファイルの最初に戻る
            print *, 'NumRec =', num_records
            if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
            do i = 1, num_records        !本読み込み
                read (Num_unit, *) matrix(:,i)
                print *, matrix(:,i)
            end do
        close (Num_unit)

    end subroutine csv_reader_int

    integer function get_records(Num_unit, header_flag)
        integer, intent(in) :: Num_unit
        logical, intent(in) :: header_flag
        integer ios

        get_records = 0
        if(header_flag) read (Num_unit, '()') !ヘッダーの読み飛ばし
        do        !レコード数を調べるループ
            read(Num_unit, '()', iostat=ios)
            if(ios/=0) exit
            get_records = get_records + 1
        end do

    end function get_records

end module csv_reader