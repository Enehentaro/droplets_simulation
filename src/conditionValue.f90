module conditionValue_m
    implicit none
    private

    type, public :: conditionValue_t
        double precision dt, L, U, direction_g(3)
        character(:), allocatable :: initialDistributionFName
        integer restart, stepEnd, outputInterval, num_drop
        real T, RH

        character(:), allocatable :: path2FlowFile
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        contains

        procedure :: read => read_condition

    end type

    contains

    subroutine read_condition(self, dir)
        use filename_mod
        class(conditionValue_t) self
        character(*), intent(in) :: dir
        integer n_unit

        double precision delta_t, L_represent, U_represent, direction_g(3)
        character(255) initialDistributionFName
        integer num_restart, n_end, outputInterval, num_droplets
        real temperature, relativeHumidity

        character(255) PATH2FlowFile
        double precision DT_FLOW
        integer OFFSET, INTERVAL_FLOW, LoopHead, LoopTail

        namelist /dropletSetting/ num_restart, n_end, delta_t, outputInterval, temperature, relativeHumidity,&
            num_droplets, direction_g, initialDistributionFName
        namelist /flowFieldSetting/ PATH2FlowFile, DT_FLOW, OFFSET, INTERVAL_FLOW, LoopHead, LoopTail, L_represent, U_represent

        initialDistributionFName = IniDistributionFName

        OPEN(newunit=n_unit, FILE=dir//'/'//conditionFName, STATUS='OLD')
            read(n_unit, nml=dropletSetting)
            read(n_unit, nml=flowFieldSetting)
        CLOSE(n_unit)

        self%dt = delta_t
        self%L = L_represent
        self%U = U_represent
        self%direction_g(:) = direction_g(:)
        self%initialDistributionFName = trim(initialDistributionFName)
        self%restart = num_restart
        self%stepEnd = n_end
        self%outputInterval = outputInterval
        self%num_drop = num_droplets
        self%T = temperature
        self%RH = relativeHumidity

        self%PATH2FlowFile = trim(PATH2FlowFile)
        self%DT_FLOW = DT_FLOW
        self%OFFSET = OFFSET
        self%INTERVAL_FLOW = INTERVAL_FLOW
        self%LoopHead = LoopHead
        self%LoopTail = LoopTail

    end subroutine

end module conditionValue_m