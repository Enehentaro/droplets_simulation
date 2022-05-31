module timeKeeper_m
    implicit none
    private

    type DateAndTime
        character(8) date
        character(10) time
        contains
        procedure :: string => DateAndTime2string
    end type

    type, public :: TimeKeeper
        private
        real startCPUtime, lastLapCPUtime
        type(DateAndTime) startDAT
        contains
        procedure startDateAndTime, erapsedTime, lapTime
    end type

    public TimeKeeper_
    public nowDateAndTime, second2HMS

    contains

    type(TimeKeeper) function TimeKeeper_()
        call cpu_time(TimeKeeper_%startCPUtime)
        TimeKeeper_%lastLapCPUtime = TimeKeeper_%startCPUtime
        call date_and_time(TimeKeeper_%startDAT%date, TimeKeeper_%startDAT%time)
    end function

    real function lapTime(self)
        class(TimeKeeper) self
        real nowCPUtime

        call cpu_time(nowCPUtime)
        lapTime = nowCPUtime - self%lastLapCPUtime

        self%lastLapCPUtime = nowCPUtime

    end function

    function startDateAndTime(self) result(str)
        class(TimeKeeper), intent(in) :: self
        character(:), allocatable :: str

        str = self%startDAT%string()

    end function

    real function erapsedTime(self)
        class(TimeKeeper), intent(in) :: self
        real nowCPUtime

        call cpu_time(nowCPUtime)
        erapsedTime = nowCPUtime - self%startCPUtime

    end function

    function nowDateAndTime() result(str)
        type(DateAndTime) dat
        character(:), allocatable :: str

        call date_and_time(dat%date, dat%time)
        str = dat%string()
        
    end function

    function DateAndTime2string(self) result(str)
        class(DateAndTime), intent(in) :: self
        character(:), allocatable :: str

        str = self%date(1:4)//'/'//self%date(5:6)//'/'//self%date(7:8)//' ' &
              //self%time(1:2)//':'//self%time(3:4)//':'//self%time(5:6)

    end function

    function second2HMS(seconds) result(str)
        real, intent(in) :: seconds
        integer h,m,s
        character(9) str

        s = int(seconds)
        h = s/3600
        s = s - h*3600
        m = s/60
        s = s - m*60
        write(str, '(i3,":",i2,":",i2)') h, m, s

    end function

end module timeKeeper_m
