module dropletMotionSimulation
    use dropletGroup_m
    implicit none

    integer outputInterval

    integer, private :: num_restart
    integer n_start, n_end

    character(:), allocatable, private :: case_dir

    type(dropletGroup) mainDroplet

    contains

    subroutine firstSet_mainDroplet
        use caseNameList_m
        integer num_initialDroplet
        character(99) fname

        case_dir = get_caseName()

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

    subroutine update_flow_check

        if(INTERVAL_FLOW <= 0) return

        call set_STEPinFLOW(Time_onSimulation(timeStep))

        if(STEPinFLOW >= NextUpdate) call update_FlowField(first=.false.)   !流れ場の更新

    end subroutine

    subroutine update_FlowField(first)
        logical, intent(in) :: first

        if(INTERVAL_FLOW <= 0) then
            call read_steadyFlowData
        else
            if(first) call set_STEPinFLOW(Time_onSimulation(timeStep))
            call read_unsteadyFlowData

        end if

        call set_MinMaxCDN

        if(first) call preprocess_onFlowField         !流れ場の前処理

        if(unstructuredGrid) then
            call boundary_setting(first)
            if(.not.first) call mainDroplet%boundary_move()
        end if
            
    end subroutine

    subroutine output_initialDroplet

        if(num_restart <= 0) call output_mainDroplet(initial = .true.)
        
    end subroutine

    subroutine output_mainDroplet(initial)
        logical, intent(in) :: initial
        character(99) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call mainDroplet%output_VTK(fname, initial)

        fname = case_dir//'/particle.csv'
        call mainDroplet%output_CSV(fname, Time_onSimulation(timeStep, dimension=.true.), initial)

        if(.not.initial) call mainDroplet%output_backup(case_dir//'/backup')

    end subroutine

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