module dropletMotionSimulation
    use dropletGroup_m
    implicit none

    integer outputInterval

    integer, private :: num_restart
    integer n_start, n_end

    character(:), allocatable :: start_date
    real start_time

    character(:), allocatable, private :: case_dir

    type(dropletGroup) mainDroplet

    contains

    subroutine firstSet_mainDroplet
        integer num_initialDroplet
        character(99) fname

        call read_and_set_condition(case_dir, num_droplet=num_initialDroplet)

        call update_FlowField(first=.true.)                !流れ場の取得

        if(num_restart <= 0) then

            if(num_restart==0) then
                call random_set  !実行時刻に応じた乱数シード設定
                mainDroplet = generate_dropletGroup(num_initialDroplet)

            else if(num_restart==-1) then
                mainDroplet = read_InitialDistribution(case_dir)

            end if

            n_start = 0
            timeStep = n_start

        else

            print*, '**RESTART**'

            write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') num_restart
            mainDroplet = read_backup(fname)   !ここで自動割り付け

            n_start = num_restart
            timeStep = n_start
            
        end if

        print*, 'num_droplets =', size(mainDroplet%droplet)

    end subroutine

    subroutine read_and_set_condition(dir, num_droplet)
        use filename_mod
        use path_operator_m
        character(*), intent(in) ::dir
        double precision dt, L, U
        double precision :: direction_g(3)
        character(99) path2FlowFile
        integer i, n_unit
        integer, optional, intent(out) :: num_droplet
        real T, RH

        OPEN(newunit=n_unit, FILE=dir//'/'//conditionFName, STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,*) num_restart
            read(n_unit,'()')
            read(n_unit,*) n_end
            read(n_unit,'()')
            read(n_unit,*) dt
            read(n_unit,'()')
            read(n_unit,*) outputInterval
            read(n_unit,'()')
            read(n_unit,*) T
            read(n_unit,*) RH
            read(n_unit,'()')
            if(present(num_droplet)) then
                read(n_unit,*) num_droplet
            else
                read(n_unit, '()')
            end if
            read(n_unit,'()')
            read(n_unit,*) (direction_g(i), i=1,3)
            
            read(n_unit,'()')

            read(n_unit,'()')
            read(n_unit,'(A)') path2FlowFile
            read(n_unit,'()')
            read(n_unit,*) DT_FLOW
            read(n_unit,'()')
            read(n_unit,*) OFFSET
            read(n_unit,'()')
            read(n_unit,*) INTERVAL_FLOW
            read(n_unit,'()')
            read(n_unit,*) LoopS
            read(n_unit,*) LoopF
            read(n_unit,'()')
            read(n_unit,*) L
            read(n_unit,*) U

        CLOSE(n_unit)

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
        print*, 'Delta_Time =', dt
        print*, 'Delta_Time inFLOW =', DT_FLOW

        call set_dir_from_path(path2FlowFile, PATH_FlowDIR, FNAME_FMT)

        call check_FILE_GRID    !気流ファイルのタイプをチェック

        call set_basical_variables(dt, L, U)

        call set_gravity_acceleration(direction_g)

        call set_environment(T, RH)

    end subroutine

    subroutine check_FlowFieldUpdate

        if(INTERVAL_FLOW <= 0) return

        call set_STEPinFLOW(TimeOnSimu())

        if(STEPinFLOW >= NextUpdate) call update_FlowField(first=.false.)   !流れ場の更新

    end subroutine

    subroutine update_FlowField(first)
        logical, intent(in) :: first

        if(INTERVAL_FLOW <= 0) then
            call read_steadyFlowData
        else
            if(first) call set_STEPinFLOW(TimeOnSimu())
            call read_unsteadyFlowData

        end if

        call set_MinMaxCDN

        if(first) call preprocess_onFlowField         !流れ場の前処理

        if(unstructuredGrid) then
            call boundary_setting(first)
            if(.not.first) call mainDroplet%boundary_move()
        end if
            
    end subroutine

    subroutine mainDroplet_process

        call mainDroplet%adhesion_check()

        call mainDroplet%survival_check()           !生存率に関する処理

        call coalescence_process        !飛沫間の合体判定

        call mainDroplet%calculation()     !飛沫の運動計算

    end subroutine

    subroutine output_mainDroplet(initial)
        logical, intent(in) :: initial
        character(99) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call mainDroplet%output_VTK(fname, initial)

        fname = case_dir//'/particle.csv'
        call mainDroplet%output_CSV(fname, TimeOnSimu(dimension=.true.), initial)

        if(.not.initial) call mainDroplet%output_backup(case_dir//'/backup')

    end subroutine

    subroutine coalescence_process
        use terminalControler_m
        integer floatings, num_coalescence
        integer, save :: last_coalescence = 0, last_floatings = 0

        floatings = mainDroplet%counter('floating')
        if(floatings > last_floatings) last_coalescence = timeStep    !浮遊数が増加したら付着判定再起動のため更新
        last_floatings = floatings

        !最後の合体から100ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((timeStep - last_coalescence) > 100) return

        call set_formatTC('(" Coalescence_check [step:", i10, "/", i10, "]")')
        call print_sameLine([timeStep, last_coalescence+100])

        call mainDroplet%coalescence_check(stat = num_coalescence)

        if(num_coalescence >= 1) last_coalescence = timeStep

    end subroutine

    subroutine checkpoint
        character(1) input
        character d_start*8, t_start*10

        if(num_restart <= 0) call output_mainDroplet(initial = .true.)

        do
            print*, 'Do you want to start the calculation? (y/n)'
            read(5,*) input

            select case(input)
                case('y')
                    exit

                case('n')
                    stop

            end select

        end do

        call cpu_time(start_time)
        call date_and_time(date = d_start, time = t_start)
        start_date = '[Start Date] ' // DateAndTime_string(d_start, t_start)

    end subroutine

    subroutine create_CaseDirectory
        use caseNameList_m
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
            write(n_unit,'(A18, F18.2)') 'Temp [degC] =', environment('Temperature')
            write(n_unit,'(A18, F18.2)') 'RH [%] =', environment('RelativeHumidity')
            write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', trim(PATH_FlowDIR)//trim(FNAME_FMT)
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