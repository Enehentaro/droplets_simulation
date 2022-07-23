module conditionValue_m
    implicit none
    private

    !>条件値クラス
    !>条件指定ファイル（namelist）を読み込んだ結果を格納する
    type, public :: conditionValue_t
        double precision dt, L, U, direction_g(3)
        character(:), allocatable :: initialDistributionFName
        integer restart, stepEnd, outputInterval, num_drop, periodicGeneration
        real T, RH

        character(:), allocatable :: path2FlowFile, meshFile
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        contains

        procedure :: read => read_condition

    end type

    contains

    subroutine read_condition(self, dir)
        use filename_m, only : conditionFName
        class(conditionValue_t) self
        character(*), intent(in) :: dir
        integer n_unit
        character(4), parameter :: None = 'None'

        double precision delta_t, L_represent, U_represent, direction_g(3)
        character(255) :: initialDistributionFName = None
        integer num_restart, n_end, outputInterval, num_droplets
        real temperature, relativeHumidity
        integer :: periodicGeneration = 0

        character(255) PATH2FlowFile
        character(255) :: meshFile = None
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        namelist /dropletSetting/ num_restart, n_end, delta_t, outputInterval, temperature, relativeHumidity,&
                                    num_droplets, direction_g, initialDistributionFName, periodicGeneration
        namelist /flowFieldSetting/ PATH2FlowFile, meshFile, DT_FLOW, OFFSET, INTERVAL_FLOW, LoopHead, LoopTail,&
                                    L_represent, U_represent

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

    end subroutine

end module conditionValue_m