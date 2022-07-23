module flow_field_m
    use unstructuredGrid_m
    implicit none
    private
    type, public, extends(FlowFieldUnstructuredGrid) :: FlowField
        private
        integer INTERVAL                             !気流データ出力間隔
        integer LoopHead, LoopTail, OFFSET
        double precision DT
        integer STEP, NextUpdate

        character(:), allocatable :: FullFileName
        character(:), allocatable :: FileNameFormat

        contains
        private

        procedure, public :: update => update_FlowField
        procedure, public ::  get_defaultFlowFileName, isUpdateTiming
        procedure, public ::  set_time => set_timeSTEPinFLOW

        procedure set_FileNameFormat, calc_NextUpdate, get_FileNumber, clamp_STEP
        procedure :: get_requiredFileName => get_requiredFlowFieldFileName

    end type

    public FlowField_

    contains

    type(FlowField) function FlowField_(&
        time, PATH2FlowFile, DeltaT, OFFSET, outputINTERVAL, LoopHead, LoopTail, meshFile)
        
        double precision, intent(in) :: time, DeltaT
        integer, intent(in) :: OFFSET, outputINTERVAL, LoopHead, LoopTail
        character(*), intent(in) :: PATH2FlowFile
        character(*), intent(in), optional :: meshFile

        call FlowField_%set_FileNameFormat(PATH2FlowFile)

        FlowField_%INTERVAL  = outputINTERVAL

        if(FlowField_%INTERVAL   <= 0) then
            print*, 'AirFlow is Steady'

        else
            print*, 'Interval of AirFlow =', FlowField_%INTERVAL 

            FlowField_%OFFSET = OFFSET
            print*, 'OFFSET =', FlowField_%OFFSET

            FlowField_%LoopHead = LoopHead
            FlowField_%LoopTail = LoopTail
            if(FlowField_%LoopTail - FlowField_%LoopHead > 0) then
                print*, 'Loop is from', FlowField_%LoopHead, 'to', FlowField_%LoopTail
            elseif(FlowField_%LoopTail - FlowField_%LoopHead == 0) then 
                print*, 'After', FlowField_%LoopTail, ', Checkout SteadyFlow'
            end if

            FlowField_%DT = DeltaT
            print*, 'Delta_Time inFLOW =', FlowField_%DT

            call FlowField_%set_time(time)

        end if

        if(present(meshFile)) then
            FlowField_%FlowFieldUnstructuredGrid = FlowFieldUnstructuredGrid_(FlowField_%get_requiredFileName(), meshFile)
        else
            FlowField_%FlowFieldUnstructuredGrid = FlowFieldUnstructuredGrid_(FlowField_%get_requiredFileName())
        end if

        call FlowField_%calc_NextUpdate()

    end function

    subroutine set_FileNameFormat(self, PATH2FlowFile)
        use path_operator_m
        class(FlowField) self
        character(*), intent(in) :: PATH2FlowFile
        character(:), allocatable :: PATH2FlowDir, prefix, suffix, FileName
        integer i_integerPart, i_, i_dot
        integer num_digit !ファイル名の整数部桁数

        self%FullFileName = PATH2FlowFile

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

            self%FileNameFormat = '("'//PATH2FlowDir//prefix//'",'//trim(digitsFormat)//',"'//suffix//'")'

            print*, 'self%FileNameFormat : ', self%FileNameFormat

        end block

    end subroutine

    function get_requiredFlowFieldFileName(self) result(FileName)
        class(FlowField) self
        character(:), allocatable :: FileName

        if(self%INTERVAL   <= 0) then
            FileName = self%FullFileName

        else
            block
                character(255) str

                write(str, self%FileNameFormat) self%get_FileNumber()
                FileName = trim(str)

            end block

        end if

    end function

    subroutine update_FlowField(self)
        class(FlowField) self

        call self%updateWithFlowFieldFile(self%get_requiredFileName())

        call self%calc_NextUpdate()

        print*, "Update_FlowField : FIN"
            
    end subroutine

    logical function isUpdateTiming(self)
        class(FlowField) self

        if(self%STEP >= self%NextUpdate .and. self%INTERVAL > 0) then
            isUpdateTiming = .true.
        else
            isUpdateTiming = .false.
        end if

    end function

    subroutine set_timeSTEPinFLOW(self, time)
        class(FlowField) self
        DOUBLE PRECISION, intent(in) :: time

        self%STEP = int(time/self%DT) + self%OFFSET   !気流計算における時刻ステップ数に相当

    end subroutine

    subroutine calc_NextUpdate(self)
        class(FlowField) self
        integer i

        i = 0
        do while(i*self%INTERVAL + self%OFFSET <= self%STEP)
            i = i + 1
        end do
        self%NextUpdate = i*self%INTERVAL + self%OFFSET

    end subroutine

    integer function get_FileNumber(self)
        class(FlowField) self

        get_FileNumber = self%OFFSET
        do while(get_FileNumber < self%STEP)
            get_FileNumber = get_FileNumber + self%INTERVAL
        end do

        get_FileNumber = get_FileNumber + self%INTERVAL     !前進評価

        call self%clamp_STEP(get_FileNumber)

    end function

    subroutine clamp_STEP(self, STEP)
        class(FlowField) self
        integer, intent(inout) :: STEP
        integer Lamda, Delta
        
        if(STEP >= self%LoopTail) then
            Lamda = self%LoopTail - self%LoopHead

            if(Lamda > 0) then
                Delta = mod(STEP - self%LoopHead, Lamda)
                STEP = self%LoopHead + Delta

            else if(Lamda == 0) then
                STEP = self%LoopTail
                self%INTERVAL = -1
                print*, '**Checkout SteadyFlow**'
            end if

        end if

    end subroutine clamp_STEP

    function get_defaultFlowFileName(self) result(fname)
        class(FlowField) self
        character(:), allocatable :: fname

        fname = self%FullFileName

    end function
    
end module flow_field_m