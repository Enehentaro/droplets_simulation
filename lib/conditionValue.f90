module conditionValue_m
    !!飛沫計算用の諸条件を取り扱う
    implicit none
    private

    type, public :: conditionValue_t
        !!飛沫計算用の諸条件をまとめた構造体
        double precision dt, L, U, direction_g(3)
        character(:), allocatable :: initialDistributionFName
        integer restart, stepEnd, outputInterval, num_drop, periodicGeneration
        real T, RH

        character(:), allocatable :: path2FlowFile, meshFile
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        contains

        procedure isInitialDistributionSpecified
            !!飛沫初期分布ファイルが指定されたか否かを返す

        procedure isMeshFileSpecified
            !!メッシュファイルが別途指定されたか否かを返す

    end type

    public read_condition

    contains

    function read_condition(dir) result(self)
        !!条件ファイルを読み込み、結果を構造体で返す。
        !!このサブルーチン実装当時、構造体をそのままnamelistにできることを知らず、わざわざ変数ひとつひとつ定義した。
        !!現在ここを変えると進行中のプロジェクト（オフィス飛沫計算など）に影響が出るおそれがあり、触れない。
        !!いつか修正したい。
        use filename_m, only : conditionFName => conditionFileName
        type(conditionValue_t) self
        character(*), intent(in) :: dir
        integer n_unit
        character(4), parameter :: None = 'None'
        double precision delta_t, L_represent, U_represent, direction_g(3)
        character(255) :: initialDistributionFName
        integer num_restart, n_end, outputInterval, num_droplets
        real temperature, relativeHumidity
        integer :: periodicGeneration = 0
        character(255) PATH2FlowFile, meshFile
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        namelist /dropletSetting/ num_restart, n_end, delta_t, outputInterval, temperature, relativeHumidity,&
                                    num_droplets, direction_g, initialDistributionFName, periodicGeneration
        namelist /flowFieldSetting/ PATH2FlowFile, meshFile, DT_FLOW, OFFSET, INTERVAL_FLOW, LoopHead, LoopTail,&
                                    L_represent, U_represent

        
        initialDistributionFName = None     !初期化
        meshFile = None                     !初期化

        OPEN(newunit=n_unit, FILE=dir//'/'//conditionFName, status='old', action='read')
            read(n_unit, nml=dropletSetting)
            read(n_unit, nml=flowFieldSetting)
        CLOSE(n_unit)

        self%dt = delta_t
        self%L = L_represent
        self%U = U_represent
        self%direction_g(:) = direction_g(:)
        if(initialDistributionFName /= None) self%initialDistributionFName = trim(initialDistributionFName)
        self%restart = num_restart
        self%stepEnd = n_end
        self%outputInterval = outputInterval
        self%num_drop = num_droplets
        self%T = temperature
        self%RH = relativeHumidity
        self%periodicGeneration = periodicGeneration

        self%PATH2FlowFile = trim(PATH2FlowFile)
        if(meshFile /= None) self%meshFile = trim(meshFile)
        self%DT_FLOW = DT_FLOW
        self%OFFSET = OFFSET
        self%INTERVAL_FLOW = INTERVAL_FLOW
        self%LoopHead = LoopHead
        self%LoopTail = LoopTail

    end function

    logical function isInitialDistributionSpecified(self)
        class(conditionValue_t), intent(in) :: self

        isInitialDistributionSpecified = allocated(self%initialDistributionFName)

    end function

    logical function isMeshFileSpecified(self)
        class(conditionValue_t), intent(in) :: self

        isMeshFileSpecified = allocated(self%meshFile)

    end function

end module conditionValue_m