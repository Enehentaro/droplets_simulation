module dropletMotionSimulation
    use dropletGroup_m
    use flow_field
    implicit none

    integer outputInterval

    integer, private :: num_restart
    integer n_start, n_end

    character(:), allocatable :: start_date
    real start_time

    character(:), allocatable, private :: case_dir

    logical, private :: startFlag = .false., adhesionSwitch = .true.

    integer, private :: last_coalescenceStep, last_numFloating, coalescenceLimit=100

    type(dropletGroup) mainDroplet

    contains

    subroutine firstSet_mainDroplet
        integer num_initialDroplet
        character(:), allocatable :: iniDisFName

        call read_and_set_condition(case_dir, num_droplet=num_initialDroplet, initialDistributionFileName=iniDisFName)
        call read_basicSetting

        call set_dropletPlacementInformation(case_dir)

        timeStep = max(num_restart, 0)                !流れ場の取得の前に必ず時刻セット

        call create_FlowField                !流れ場の取得

        if(num_restart <= 0) then

            if(num_restart==0) then
                call random_set  !実行時刻に応じた乱数シード設定
                mainDroplet = generate_dropletGroup(num_initialDroplet, outputDir=case_dir)

            else if(num_restart==-1) then
                mainDroplet = read_InitialDistribution(case_dir//'/'//iniDisFName)

            end if

            n_start = 0

        else

            print*, '**RESTART**'

            block
                character(255) fname
                write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') num_restart
                mainDroplet = read_backup(trim(fname))   !ここで自動割り付け
            end block

            n_start = num_restart

        end if

        print*, 'num_droplets =', size(mainDroplet%droplet)
        last_coalescenceStep = 0
        last_numFloating = 0

    end subroutine

    subroutine read_and_set_condition(dir, num_droplet, initialDistributionFileName)
        use dropletEquation_m
        use filename_mod
        character(*), intent(in) :: dir
        double precision delta_t, L_represent, U_represent
        double precision :: direction_g(3)
        character(255) path2FlowFile, initialDistributionFName
        character(:), allocatable, optional, intent(out) :: initialDistributionFileName
        integer n_unit, num_droplets
        integer, optional, intent(out) :: num_droplet
        real temperature, relativeHumidity
        namelist /dropletSetting/ num_restart, n_end, delta_t, outputInterval, temperature, relativeHumidity,&
            num_droplets, direction_g, initialDistributionFName
        namelist /flowFieldSetting/ path2FlowFile, DT_FLOW, OFFSET, INTERVAL_FLOW, LoopS, LoopF, L_represent, U_represent

        initialDistributionFName = IniDistributionFName

        OPEN(newunit=n_unit, FILE=dir//'/'//conditionFName, STATUS='OLD')
            read(n_unit, nml=dropletSetting)
            read(n_unit, nml=flowFieldSetting)
        CLOSE(n_unit)

        if(present(num_droplet)) num_droplet = num_droplets
        if(present(initialDistributionFileName)) initialDistributionFileName = trim(initialDistributionFName)

        if(num_restart > 0) then
            print*, 'Restart from', num_restart
        else if(num_restart == -1) then
            print*, 'InitialDistributon is Specified.'
        end if

        print*, 'n_end =',n_end
        print*, 'output interval =', outputInterval

        if(INTERVAL_FLOW > 0) then
            print*, 'Interval of AirFlow =', INTERVAL_FLOW
        else
            print*, 'AirFlow is Steady'
        end if

        if(loopf - loops > 0) then
            print*, 'Loop is from', loops, 'to', loopf
        elseif(loopf - loops == 0) then 
            print*, 'After', loopf, ', Checkout SteadyFlow'
        end if
        print*, 'Delta_Time =', delta_t
        print*, 'Delta_Time inFLOW =', DT_FLOW

        call set_FlowFileNameFormat(path2FlowFile)

        call set_basical_variables(delta_t, L_represent, U_represent)

        call set_gravity_acceleration(direction_g)

        call set_dropletEnvironment(temperature, relativeHumidity)

    end subroutine

    subroutine read_basicSetting
        integer n_unit
        namelist /basicSetting/ coalescenceLimit, adhesionSwitch

        open(newunit=n_unit, file='option/basicSetting.nml', status='old')
            read(n_unit, nml=basicSetting)
        close(n_unit)

    end subroutine

    subroutine check_FlowFieldUpdate

        if((INTERVAL_FLOW <= 0).or.(timeStep==n_end)) return

        call set_STEPinFLOW(TimeOnSimu())

        if(isUpdateTiming()) call update_FlowField   !流れ場の更新

    end subroutine

    subroutine create_FlowField

        if(INTERVAL_FLOW <= 0) then
            call read_steadyFlowData
        else
            call set_STEPinFLOW(TimeOnSimu())
            call read_unsteadyFlowData

        end if

        call set_MinMaxCDN

        call preprocess_onFlowField         !流れ場の前処理
            
    end subroutine

    subroutine update_FlowField

        call read_unsteadyFlowData

        ! if(unstructuredGrid) then
            call boundary_setting(first=.false.)
            call mainDroplet%boundary_move()
        ! end if

        call set_MinMaxCDN
            
    end subroutine

    subroutine mainDroplet_process

        if(adhesionSwitch) call mainDroplet%adhesion_check()

        call mainDroplet%survival_check()           !生存率に関する処理

        call coalescence_process        !飛沫間の合体判定

        call mainDroplet%calculation()     !飛沫の運動計算

    end subroutine

    subroutine output_mainDroplet(initial)
        logical, intent(in) :: initial
        character(255) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call mainDroplet%output_VTK(fname, deadline=initial)

        fname = case_dir//'/particle.csv'
        call mainDroplet%output_CSV(fname, TimeOnSimu(dimension=.true.), initial)

        call mainDroplet%output_backup(case_dir//'/backup', initial)

    end subroutine

    subroutine coalescence_process
        use terminalControler_m
        integer numFloating, num_coalescence
        
        numFloating = mainDroplet%counter('floating')
        if(numFloating > last_numFloating) last_coalescenceStep = timeStep    !浮遊数が増加したら付着判定再起動のため更新
        last_numFloating = numFloating

        !最後の合体から100ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((timeStep - last_coalescenceStep) > coalescenceLimit) return

        call set_formatTC('(" Coalescence_check [step:", i10, "/", i10, "]")')
        call print_sameLine([timeStep, last_coalescenceStep + coalescenceLimit])

        call mainDroplet%coalescence_check(stat = num_coalescence)

        if(num_coalescence >= 1) last_coalescenceStep = timeStep

    end subroutine

    subroutine checkpoint
        character(1) input
        character d_start*8, t_start*10

        if(num_restart <= 0) call output_mainDroplet(initial = .true.)

        if(.not.startFlag) then
            do
                print*, 'Do you want to start the calculation? (y/n)'
                read(5,*) input

                select case(input)
                    case('y')
                        startFlag = .true.
                        exit

                    case('n')
                        stop

                end select

            end do
        end if

        call cpu_time(start_time)
        call date_and_time(date = d_start, time = t_start)
        start_date = '[Start Date] ' // DateAndTime_string(d_start, t_start)

    end subroutine

    subroutine create_CaseDirectory
        use caseName_m
        use path_operator_m

        case_dir = get_caseName()

        print*, '#', nowCase, '[',case_dir,']'

        call make_directory(case_dir//'/VTK')
        call make_directory(case_dir//'/backup')
        
    end subroutine

    subroutine periodicOutput
        use terminalControler_m

        print*, start_date
        print*, 'Now_Step_Time =', TimeOnSimu(dimension=.true.), '[sec]'
        print*, '# floating :', mainDroplet%Counter('floating'), '/', mainDroplet%Counter('total')
        if(refCellSearchInfo('FalseRate') >= 1) print*, '# searchFalse :', refCellSearchInfo('NumFalse')
        call output_mainDroplet(initial=.false.)
        print '("====================================================")'
        call reset_formatTC

    end subroutine

    subroutine output_ResultSummary()
        use dropletEquation_m
        integer n_unit, cnt
        real end_time
        character(50) fname
        character d_end*8, t_end*10
        logical existance
        double precision TimeStart, TimeEnd
        character(:), allocatable :: end_date
        
        call cpu_time(end_time)
        call date_and_time(date = d_end, time = t_end)

        end_date = '[ END  Date] ' // DateAndTime_string(d_end, t_end)
        print*, start_date
        print*, end_date

        fname = case_dir//'/ResultSummary.txt'
        inquire(file=fname, exist=existance)
        cnt = 0
        do while(existance)
            cnt = cnt + 1
            write(fname,'("'//case_dir//'/ResultSummary_", i0, ".txt")') cnt
            inquire(file=fname, exist=existance)
        end do

        TimeStart = TimeOnSimu(step=n_start, dimension=.true.)
        TimeEnd = TimeOnSimu(step=n_end, dimension=.true.)

        open(newunit=n_unit, file=fname, status='new')
            write(n_unit,*)'*******************************************'
            write(n_unit,*)'*                                         *'
            write(n_unit,*)'*             Result Summary              *'
            write(n_unit,*)'*                                         *'
            write(n_unit,*)'*******************************************'
            write(n_unit,'(A)') '======================================================='
            write(n_unit,*) start_date
            write(n_unit,*) end_date
            write(n_unit, '(A18, F15.3, 2X, A)') 'Erapsed Time =', end_time - start_time, '[sec]'
            write(n_unit, '(A18, F15.3, 2X, A)') 'Cost of Calc =', &
                    (end_time - start_time) / (TimeEnd - TimeStart), '[sec/sec]'
            write(n_unit,'(A)') '======================================================='
            write(n_unit, '(A18, 2(F15.3,2x,A))') 'Time [sec] =', TimeStart, '-', TimeEnd
            write(n_unit, '(A18, 2(I15,2x,A))') 'Step =', n_start, '-', n_end !計算回数
            write(n_unit,'(A18, I15)') 'OutputInterval =', outputInterval
            write(n_unit,'(A)') '======================================================='
            write(n_unit,'(A18, I15)') '#Droplets =', mainDroplet%counter('total')
            write(n_unit,'(A18, I15)') 'floating =', mainDroplet%counter('floating')
            write(n_unit,'(A18, I15)') 'death =', mainDroplet%counter('death') !生存率で消滅
            write(n_unit,'(A18, I15)') 'coalescence =', mainDroplet%counter('coalescence') !生存率で消滅
            write(n_unit,'(A18, I15)') 'adhesion =', mainDroplet%counter('adhesion') !付着したすべてのウイルス数
            write(n_unit,'(A)') '======================================================='
            write(n_unit,'(A18, F18.2)') 'Temp [degC] =', dropletEnvironment('Temperature')
            write(n_unit,'(A18, F18.2)') 'RH [%] =', dropletEnvironment('RelativeHumidity')
            write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', get_FlowFileName()
            write(n_unit, '(A18, 2(I15,2x,A))') 'SearchFalseInfo :', refCellSearchInfo('NumFalse'), &
                    ' (', refCellSearchInfo('FalseRate'), '%)'

        close(n_unit)
        
    end subroutine

    function DateAndTime_string(date, time) result(string)
        character(*), intent(in) :: date, time
        character(:), allocatable :: string

        string = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' ' &
              //time(1:2)//':'//time(3:4)//':'//time(5:6)

    end function

    subroutine random_set    !実行時刻に依存した乱数シードを指定する
        implicit none
        integer :: seedsize, i
        integer, allocatable :: seed(:)

        print*, 'call:random_set'
    
        call random_seed(size=seedsize) !シードのサイズを取得。（コンパイラごとに異なるらしい）
        allocate(seed(seedsize)) !新シード配列サイズの割り当て
    
        do i = 1, seedsize
            call system_clock(count=seed(i)) !時間を新シード配列に取得
        end do
    
        call random_seed(put=seed(:)) !新シードを指定
          
    end subroutine random_set

end module dropletMotionSimulation