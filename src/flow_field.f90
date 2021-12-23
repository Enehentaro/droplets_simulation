module flow_field
    use unstructuredGrid_mod
    ! use CUBE_mod
    implicit none

    integer INTERVAL_FLOW                           !気流データ出力間隔
    integer LoopHead, LoopTail, OFFSET
    double precision DT_FLOW
    integer, private :: STEPinFLOW, NextUpdate

    character(:), allocatable, private :: PATH_FlowDIR, HEAD_AIR, FNAME_FMT !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    integer, private :: FNAME_DIGITS !ファイル名の整数部桁数
    
    logical, private :: unstructuredGrid

    real, private :: MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限

    ! type reference_cell_t
    !     integer :: ID = 0, nodeID(3) = 0
    ! end type

    integer, private :: num_refCellSearchFalse, num_refCellSearch

    contains

    subroutine check_FILE_GRID
        integer i_integerPart, i_

        unstructuredGrid = .true.

        i_integerPart = index(FNAME_FMT, '0')   !ひとまず最初のゼロの位置を整数部位置とする
        i_ = index(FNAME_FMT, '_', back=.true.)  !アンダーバーの位置取得
        if(i_ > i_integerPart) i_integerPart = i_ + 1   !アンダーバーの位置がゼロの位置より後ろの場合、アンダーバー以降を整数部位置とする
        HEAD_AIR = FNAME_FMT(: i_integerPart - 1)             !ファイル名の接頭部(整数部位置の手前まで)
        FNAME_DIGITS = index(FNAME_FMT, '.') - i_integerPart   !ファイル名の整数部桁数(整数部位置からドットまでの文字数)

        if(FNAME_FMT(len(FNAME_FMT)-1 :) == '.f') then
            unstructuredGrid = .false.
            print*, 'FILE_GRID : [CUBE] ', trim(FNAME_FMT)
            print*, 'DO NOT USE CUBEGRID'
            stop

        else
            if(index(FNAME_FMT, ',') > 0) then  !カンマが存在する場合
                call check_FILE_TYPE(FNAME_FMT(:i_-1), PATH_FlowDIR//FNAME_FMT(i_+1:))
                FNAME_FMT = FNAME_FMT(1:i_-1)
            else
                call check_FILE_TYPE(FNAME_FMT)
            end if

        end if

    end subroutine check_FILE_GRID

    function get_digits_format() result(format)
        character(2*(FNAME_DIGITS/10 +1) + 2) format

        if(FNAME_DIGITS <= 9) then
            write(format,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS
        else
            write(format,'("i", i2, ".", i2)') FNAME_DIGITS, FNAME_DIGITS
        end if

    end function get_digits_format

    subroutine preprocess_onFlowField
        ! use adhesion_onSTL_m
        use adjacent_information
        logical success

        ! if(unstructuredGrid) then
            call read_adjacency(PATH_FlowDIR, success)
            if(success) then
                call read_boundaries(PATH_FlowDIR)

            else
                call solve_adjacentInformation
                call output_boundaries(PATH_FlowDIR)
                call output_adjacency(PATH_FlowDIR)

            end if

            call boundary_setting(first=.true.)

        ! else
        !     call read_faceShape(PATH_FlowDIR)
        !     call set_faceShape
        ! end if

        num_refCellSearchFalse = 0
        num_refCellSearch = 0

    end subroutine preprocess_onFlowField

    subroutine create_FlowField(time, PATH2FlowFile, DeltaT_FLOW, timeOFFSET, outputINTERVAL_FLOW, flowLoopHead, flowLoopTail)
        use path_operator_m
        double precision, intent(in) :: time, DeltaT_FLOW
        integer, intent(in) :: timeOFFSET, outputINTERVAL_FLOW, flowLoopHead, flowLoopTail
        character(*), intent(in) :: PATH2FlowFile

        call set_dir_from_path(PATH2FlowFile, PATH_FlowDIR, FNAME_FMT)
        call check_FILE_GRID    !気流ファイルのタイプをチェック

        INTERVAL_FLOW = outputINTERVAL_FLOW

        if(INTERVAL_FLOW <= 0) then
            print*, 'AirFlow is Steady'
            call read_steadyFlowData
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
            call read_unsteadyFlowData

        end if

        call set_MinMaxCDN

        call preprocess_onFlowField         !流れ場の前処理
            
    end subroutine

    subroutine read_steadyFlowData
        character(:), allocatable :: FNAME

        ! if(unstructuredGrid) then
            FNAME = trim(PATH_FlowDIR)//trim(FNAME_FMT)
            call read_unstructuredGrid_byNAME(FNAME)
            call set_gravity_center

        ! else
        !     FNAME = trim(PATH_FlowDIR)//trim(FNAME_FMT)
        !     call read_CUBE_data(FNAME, trim(PATH_FlowDIR))

        ! end if
            
    end subroutine

    subroutine read_unsteadyFlowData
        integer FNUM
        character(:), allocatable :: FNAME, digits_fmt

        FNUM = get_FileNumber()
        digits_fmt = get_digits_format()

        ! if(unstructuredGrid) then
            FNAME = trim(PATH_FlowDIR)//trim(HEAD_AIR)
            call read_unstructuredGrid_byNumber(FNAME, digits_fmt, FNUM)
            call set_gravity_center

        ! else
        !     write(FNAME,'("'//trim(PATH_FlowDIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".f")') FNUM
        !     call read_CUBE_data(FNAME, trim(PATH_FlowDIR))

        ! end if

        call calc_NextUpdate
            
    end subroutine

    subroutine set_MinMaxCDN
        ! real min_max(6)

        ! if(unstructuredGrid) then
            MAX_CDN(1) = maxval(NODEs(:)%coordinate(1))
            MAX_CDN(2) = maxval(NODEs(:)%coordinate(2))
            MAX_CDN(3) = maxval(NODEs(:)%coordinate(3))
            print*, 'MAX_coordinates=', MAX_CDN(:)
                
            MIN_CDN(1) = minval(NODEs(:)%coordinate(1))
            MIN_CDN(2) = minval(NODEs(:)%coordinate(2))
            MIN_CDN(3) = minval(NODEs(:)%coordinate(3))
            print*, 'MIN_coordinates=', MIN_CDN(:)

        ! else
            ! min_max = get_minMax_CUBE()

            ! MIN_CDN(:) = min_max(1:3)
            ! MAX_CDN(:) = min_max(4:6)
        ! end if      

    end subroutine set_MinMaxCDN

    function get_areaMinMax() result(MinMax)
        real MinMax(3,2)

        MinMax(:,1) = MIN_CDN
        MinMax(:,2) = MAX_CDN

    end function

    subroutine search_refCELL(X, reference_cell, stat)
        real, intent(in) :: X(3)
        integer, intent(inout) :: reference_cell
        logical, optional :: stat

        num_refCellSearch = num_refCellSearch + 1

        reference_cell = nearer_cell(X, reference_cell)
        if(present(stat)) stat = .True.

        if (.not.nearcell_check(X(:), reference_cell)) then
            reference_cell = nearest_cell(X)
            if(present(stat)) stat = .false.
            num_refCellSearchFalse = num_refCellSearchFalse + 1
        end if
    
    end subroutine search_refCELL

    ! subroutine search_refCELL_onCUBE(X, reference_cell)
    !     real, intent(in) :: X(3)
    !     type(reference_cell_t), intent(inout) :: reference_cell

    !     num_refCellSearch = num_refCellSearch + 1

    !     reference_cell%nodeID(:) = nearer_node(X, reference_cell%nodeID, reference_cell%ID)
    !     if (.not.nearNode_check(X, reference_cell%nodeID, reference_cell%ID)) then
    !         reference_cell%ID = get_cube_contains(X)    
    !         reference_cell%nodeID(:) = nearest_node(X, reference_cell%ID)
    !         num_refCellSearchFalse = num_refCellSearchFalse + 1
    !     end if

    ! end subroutine search_refCELL_onCUBE

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

        i = 0
        do while(i*INTERVAL_FLOW + OFFSET <= STEPinFLOW)
            i = i + 1
        end do

        NextUpdate = i*INTERVAL_FLOW + OFFSET

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

    function get_FlowFileName() result(fname)
        character(:), allocatable :: fname

        fname = PATH_FlowDIR//FNAME_FMT

    end function

    integer function refCellSearchInfo(name)
        character(*), intent(in) :: name

        select case(name)
            case('NumFalse')
                refCellSearchInfo = num_refCellSearchFalse
            case('FalseRate')
                refCellSearchInfo = 100 * num_refCellSearchFalse / (num_refCellSearch + 1)
            case default
                print*, 'ERROR refCellSearchInfo : ', name
                stop
        end select

    end function
    
end module flow_field