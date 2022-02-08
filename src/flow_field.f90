module flow_field
    use unstructuredGrid_mod
    implicit none

    integer INTERVAL_FLOW                           !気流データ出力間隔
    integer LoopHead, LoopTail, OFFSET
    double precision DT_FLOW
    integer, private :: STEPinFLOW, NextUpdate

    character(:), allocatable, private :: FullFileName
    character(:), allocatable, private :: FileNameFormat

    type(UnstructuredGridAdjacencySolved) mainMesh

    contains

    subroutine set_FileNameFormat(PATH2FlowFile)
        use path_operator_m
        character(*), intent(in) :: PATH2FlowFile
        character(:), allocatable :: PATH2FlowDir, prefix, suffix, FileName
        integer i_integerPart, i_, i_dot
        integer num_digit !ファイル名の整数部桁数

        FullFileName = PATH2FlowFile

        call get_DirFromPath(PATH2FlowFile, PATH2FlowDir, FileName)

        i_integerPart = index(FileName, '0')   !ひとまず最初のゼロの位置を整数部位置とする
        i_ = index(FileName, '_', back=.true.)  !アンダーバーの位置取得
        if(i_ > i_integerPart) i_integerPart = i_ + 1   !アンダーバーの位置がゼロの位置より後ろの場合、アンダーバー以降を整数部位置とする

        prefix = FileName(: i_integerPart - 1)             !ファイル名の接頭部(整数部位置の手前まで)

        i_dot = index(FileName, '.')   !ドットの位置
        num_digit = i_dot - i_integerPart   !ファイル名の整数部桁数(整数部位置からドットまでの文字数)

        suffix = FileName(i_dot : )

        block
            character(2*(num_digit/10 +1) + 2) digitsFormat

            if(suffix=='.fld') then
                digitsFormat = 'i0'
            else
                write(digitsFormat,'("i", i0, ".", i0)') num_digit, num_digit
            end if

            FileNameFormat = '("'//PATH2FlowDir//prefix//'",'//trim(digitsFormat)//',"'//suffix//'")'

            print*, 'FileNameFormat : ', FileNameFormat

        end block

    end subroutine

    function get_requiredFlowFieldFileName() result(FileName)
        character(:), allocatable :: FileName

        if(INTERVAL_FLOW <= 0) then
            FileName = FullFileName

        else
            block
                character(255) str

                write(str, FileNameFormat) get_FileNumber()
                FileName = trim(str)

            end block

        end if

    end function

    subroutine create_FlowField(time, PATH2FlowFile, DeltaT_FLOW, timeOFFSET, outputINTERVAL_FLOW, flowLoopHead, flowLoopTail, &
                                 meshFile)
        double precision, intent(in) :: time, DeltaT_FLOW
        integer, intent(in) :: timeOFFSET, outputINTERVAL_FLOW, flowLoopHead, flowLoopTail
        character(*), intent(in) :: PATH2FlowFile
        character(*), intent(in), optional :: meshFile

        call set_FileNameFormat(PATH2FlowFile)

        INTERVAL_FLOW = outputINTERVAL_FLOW

        if(INTERVAL_FLOW <= 0) then
            print*, 'AirFlow is Steady'

        else
            print*, 'Interval of AirFlow =', INTERVAL_FLOW

            OFFSET = timeOFFSET
            print*, 'OFFSET =', OFFSET

            LoopHead = flowLoopHead
            LoopTail = flowLoopTail
            if(LoopTail - LoopHead > 0) then
                print*, 'Loop is from', LoopHead, 'to', LoopTail
            elseif(LoopTail - LoopHead == 0) then 
                print*, 'After', LoopTail, ', Checkout SteadyFlow'
            end if

            DT_FLOW = DeltaT_FLOW
            print*, 'Delta_Time inFLOW =', DT_FLOW

            call set_STEPinFLOW(time)

        end if

        if(present(meshFile)) then
            mainMesh = UnstructuredGridAdjacencySolved_(get_requiredFlowFieldFileName(), meshFile)
        else
            mainMesh = UnstructuredGridAdjacencySolved_(get_requiredFlowFieldFileName())
        end if

        call calc_NextUpdate

    end subroutine

    subroutine update_FlowField

        call mainMesh%updateWithFlowFieldFile(get_requiredFlowFieldFileName())

        call calc_NextUpdate
            
    end subroutine

    logical function isUpdateTiming()

        if(STEPinFLOW >= NextUpdate) then
            isUpdateTiming = .true.
        else
            isUpdateTiming = .false.
        end if

    end function

    subroutine set_STEPinFLOW(time)
        DOUBLE PRECISION, intent(in) :: time

        STEPinFLOW = int(time/DT_FLOW) + OFFSET   !気流計算における時刻ステップ数に相当

    end subroutine

    subroutine calc_NextUpdate
        integer i

        if(INTERVAL_FLOW > 0) then
            i = 0
            do while(i*INTERVAL_FLOW + OFFSET <= STEPinFLOW)
                i = i + 1
            end do
            NextUpdate = i*INTERVAL_FLOW + OFFSET

        else
            NextUpdate = 999999999

        end if



    end subroutine

    integer function get_FileNumber()

        get_FileNumber = OFFSET
        do while(get_FileNumber < STEPinFLOW)
            get_FileNumber = get_FileNumber + INTERVAL_FLOW
        end do

        get_FileNumber = get_FileNumber + INTERVAL_FLOW   !前進評価

        call clamp_STEP(get_FileNumber)

    end function

    subroutine clamp_STEP(STEP)
        integer, intent(inout) :: STEP
        integer Lamda, Delta
        
        if(STEP >= LoopTail) then
            Lamda = LoopTail - LoopHead

            if(Lamda > 0) then
                Delta = mod(STEP - LoopHead, Lamda)
                STEP = LoopHead + Delta

            else if(Lamda == 0) then
                STEP = LoopTail
                INTERVAL_FLOW = -1
                print*, '**Checkout SteadyFlow**'
            end if
        end if

    end subroutine clamp_STEP

    function get_defaultFlowFileName() result(fname)
        character(:), allocatable :: fname

        fname = FullFileName

    end function
    
end module flow_field