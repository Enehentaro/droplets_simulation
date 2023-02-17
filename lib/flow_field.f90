module flow_field_m
    use unstructuredGrid_m
    implicit none
    private
    type, public, extends(FlowFieldUnstructuredGrid) :: FlowField
        private
        integer INTERVAL
            !! 気流ファイル出力間隔

        integer LoopHead
            !! 流れ場ファイルのループの先頭

        integer LoopTail
            !! 流れ場ファイルのループの末尾

        integer OFFSET
            !! 流れ場と飛沫の間の時間ステップオフセット

        double precision DT
            !! 流体計算の時間間隔

        integer STEP
            !! 時間ステップ

        integer NextUpdate
            !! 次に流れ場更新を行うべき時間ステップ

        character(:), allocatable :: FullFileName
            !! デフォルトの流れ場ファイル名

        character(:), allocatable :: FileNameFormat
            !! 流れ場ファイルのファイル名のフォーマット文字列

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

    function FlowField_(time, PATH2FlowFile, DeltaT, OFFSET, outputINTERVAL, LoopHead, LoopTail, meshFile)&
        result(flow_field)
        
        double precision, intent(in) :: time
            !! 現在無次元時刻

        double precision, intent(in) :: DeltaT
            !! 流体計算の時間間隔

        integer, intent(in) :: OFFSET, outputINTERVAL, LoopHead, LoopTail
        character(*), intent(in) :: PATH2FlowFile
        character(*), intent(in), optional :: meshFile
        type(FlowField) flow_field

        call flow_field%set_FileNameFormat(PATH2FlowFile)

        flow_field%INTERVAL  = outputINTERVAL

        if(flow_field%INTERVAL   <= 0) then
            print*, 'AirFlow is Steady'

        else
            print*, 'Interval of AirFlow =', flow_field%INTERVAL 

            flow_field%OFFSET = OFFSET
            print*, 'OFFSET =', flow_field%OFFSET

            flow_field%LoopHead = LoopHead
            flow_field%LoopTail = LoopTail
            if(flow_field%LoopTail - flow_field%LoopHead > 0) then
                print*, 'Loop is from', flow_field%LoopHead, 'to', flow_field%LoopTail
            elseif(flow_field%LoopTail - flow_field%LoopHead == 0) then 
                print*, 'After', flow_field%LoopTail, ', Checkout SteadyFlow'
            end if

            flow_field%DT = DeltaT
            print*, 'Delta_Time inFLOW =', flow_field%DT

            call flow_field%set_time(time)

        end if

        if(present(meshFile)) then
            flow_field%FlowFieldUnstructuredGrid = FlowFieldUnstructuredGrid_withMeshFile(flow_field%get_requiredFileName(), meshFile)
        else
            flow_field%FlowFieldUnstructuredGrid = FlowFieldUnstructuredGrid_(flow_field%get_requiredFileName())
        end if

        call flow_field%calc_NextUpdate()

    end function

    subroutine set_FileNameFormat(self, PATH2FlowFile)
        !! ファイル名から連番の位置を取得し、ファイル名のフォーマットを推論する
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
        !! 現時刻ステップに対応する流れ場ファイル名を返す
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
        !! 流れ場を更新する
        class(FlowField) self

        call self%updateWithFlowFieldFile(self%get_requiredFileName())

        call self%calc_NextUpdate()

        print*, "Update_FlowField : FIN"
            
    end subroutine

    logical function isUpdateTiming(self)
        !! 現在ステップが流れ場更新のタイミングかどうかを判定
        class(FlowField) self

        if(self%STEP >= self%NextUpdate .and. self%INTERVAL > 0) then
            !現在ステップが更新ステップ以降でなおかつ非定常流
            isUpdateTiming = .true.
        else
            isUpdateTiming = .false.
        end if

    end function

    subroutine set_timeSTEPinFLOW(self, time)
        !! 流れ場に現在時刻を通知
        class(FlowField) self
        DOUBLE PRECISION, intent(in) :: time
            !! 飛沫計算における時刻

        self%STEP = int(time/self%DT) + self%OFFSET   !気流計算における時刻ステップ数に相当

    end subroutine

    subroutine calc_NextUpdate(self)
        !! 次の更新時間ステップを算出
        class(FlowField) self
        integer i

        i = 0
        do while(i*self%INTERVAL + self%OFFSET <= self%STEP)
            i = i + 1
        end do
        self%NextUpdate = i*self%INTERVAL + self%OFFSET

    end subroutine

    function get_FileNumber(self) result(FileNumber)
        !! 現在時刻に対応する流れ場連番ファイル番号を計算
        class(FlowField) self
        integer FileNumber

        FileNumber = self%OFFSET
        do while(FileNumber < self%STEP)
            FileNumber = FileNumber + self%INTERVAL
        end do

        FileNumber = FileNumber + self%INTERVAL     !前進評価

        call self%clamp_STEP(FileNumber)

    end function

    subroutine clamp_STEP(self, STEP)
        !! 連番ファイルをループさせる場合を想定し、ファイル番号をclamp
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

    end subroutine

    function get_defaultFlowFileName(self) result(fname)
        !! デフォルトの流れ場ファイル名を返す
        class(FlowField) self
        character(:), allocatable :: fname

        fname = self%FullFileName

    end function
    
end module flow_field_m