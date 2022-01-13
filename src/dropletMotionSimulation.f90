module dropletMotionSimulation
    use dropletGenerator_m
    use dropletEquation_m
    use flow_field
    implicit none

    private

    integer, public :: outputInterval

    integer num_restart
    integer, public, target :: timeStep
    integer, public :: n_start, n_end

    character(:), allocatable :: start_date
    real start_time

    character(:), allocatable :: case_dir

    logical :: startFlag = .false., adhesionSwitch = .true.

    integer :: last_coalescenceStep, last_numFloating, coalescenceLimit=100, num_divide=4

    type(DropletGroup) mainDroplet

    type(DropletEquationSolver), target :: dropletSolver

    type(DropletGenerator) dropGenerator

    public create_CaseDirectory, firstSet_mainDroplet, checkpoint, output_ResultSummary
    public mainDroplet_process, periodicOutput, check_FlowFieldUpdate

    contains

    subroutine firstSet_mainDroplet
        use conditionValue_m
        type(conditionValue_t) condVal

        call condVal%read(case_dir)
        num_restart = condVal%restart
        n_end = condVal%stepEnd
        outputInterval = condVal%outputInterval
        print*, 'n_end =',n_end
        print*, 'output interval =', outputInterval

        dropletSolver = DropletEquationSolver_( &
                            condVal%dt, condVal%L, condVal%U, &
                            condVal%direction_g, condVal%T, condVal%RH &
                            )

        n_start = max(num_restart, 0)
        ! dropGenerator%timeStep = n_start              !流れ場の取得の前に必ず時刻セット

        call create_FlowField( &                !流れ場の取得
            dropletSolver%TimeOnSimu(n_start), condVal%PATH2FlowFile, &
            condVal%DT_FLOW, condVal%OFFSET, condVal%INTERVAL_FLOW, &
            condVal%LoopHead, condVal%LoopTail &
            )

        block
            character(:), allocatable :: radDisFNAME

            call read_basicSetting(radiusDistributionFilename = radDisFNAME)

            dropGenerator = DropletGenerator_(dropletSolver, radDisFNAME, case_dir)

        end block

        if(num_restart <= 0) then

            if(num_restart==0) then
                call random_set  !実行時刻に応じた乱数シード設定
                mainDroplet = dropGenerator%generateDroplet(condVal%num_drop, timeStep, outputDir=case_dir)

            else if(num_restart==-1) then
                mainDroplet = read_InitialDistribution(case_dir//'/'//condVal%initialDistributionFName)

            end if

        else

            print*, '**RESTART**'

            block
                use virusDroplet_m
                character(255) fname

                write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') num_restart
                mainDroplet = read_backup(trim(fname))   !ここで自動割り付け
            end block

        end if

        print*, 'num_droplets =', size(mainDroplet%droplet)
        last_coalescenceStep = 0
        last_numFloating = 0

    end subroutine

    subroutine read_basicSetting(radiusDistributionFilename)
        integer n_unit
        character(:), allocatable, intent(out) :: radiusDistributionFilename
        character(255) radiusDistributionFNAME
        namelist /basicSetting/ coalescenceLimit, adhesionSwitch, num_divide, radiusDistributionFNAME

        open(newunit=n_unit, file='option/basicSetting.nml', status='old')
            read(n_unit, nml=basicSetting)
        close(n_unit)

        radiusDistributionFilename = trim(radiusDistributionFNAME)

    end subroutine

    type(DropletGroup) function read_InitialDistribution(fname)
        use virusDroplet_m
        character(*), intent(in) :: fname

        ! read_InitialDistribution = read_droplet_VTK(dir//'/InitialDistribution.vtk')
        ! read_InitialDistribution%droplet(:)%initialRadius = read_InitialDistribution%droplet(:)%radius
        ! call read_InitialDistribution%calc_minimumRadius(RH)

        ! call read_InitialDistribution%first_refCellSearch()

        read_initialDistribution = read_backup(fname)

        call read_initialDistribution%first_refCellSearch()

    end function

    subroutine check_FlowFieldUpdate

        if((INTERVAL_FLOW <= 0).or.(timeStep==n_end)) return

        call set_STEPinFLOW(dropletSolver%TimeOnSimu(timeStep))

        if(isUpdateTiming()) call update_FlowField   !流れ場の更新

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

        call mainDroplet%survival_check(dropletSolver%TimeOnSimu(timeStep))           !生存率に関する処理

        call divideAreacoalescence_process        !飛沫間の合体判定

        call Calculation_Droplets     !飛沫の運動計算

    end subroutine

    subroutine Calculation_Droplets()
        integer vn, targetID

        ! if(unstructuredGrid) then
            !$omp parallel do
            do vn = 1, size(mainDroplet%droplet)

                select case(mainDroplet%droplet(vn)%status)
                case(0)
                    call evaporation(mainDroplet%droplet(vn))    !蒸発方程式関連の処理
                    call motionCalculation(mainDroplet%droplet(vn))     !運動方程式関連の処理

                case(-2)
                    targetID = mainDroplet%droplet(vn)%coalesID  !合体飛沫の片割れも移動させる
                    mainDroplet%droplet(vn)%position = mainDroplet%droplet(targetID)%position
                    mainDroplet%droplet(vn)%velocity = mainDroplet%droplet(targetID)%velocity

                end select

            end do
            !$omp end parallel do

        ! else
        !     do vn = 1, size(self%droplet)
        !         if(self%droplet(vn)%status /= 0) cycle !浮遊状態でないなら無視
        !         call self%droplet(vn)%evaporation()    !蒸発方程式関連の処理
        !         call self%droplet(vn)%motionCalculation_onCUBE()     !運動方程式関連の処理
        !     end do

        ! end if

    end subroutine

    subroutine evaporation(droplet) !CALCULATE droplet evaporation
        use virusDroplet_m
        type(virusDroplet_t) droplet
        double precision radius_n
      
        if (droplet%radius <= droplet%radius_min) return  !半径が最小になったものを除く
    
        radius_n = dropletSolver%evaporatin_eq(droplet%radius)
        
        droplet%radius = max(radius_n, droplet%radius_min)
      
    end subroutine


    subroutine motionCalculation(droplet)
        use virusDroplet_m
        type(virusDroplet_t) droplet
        double precision velAir(3)

        velAir(:) = CELLs(droplet%refCellID)%flowVelocity(:)
    
        call dropletSolver%solve_motionEquation(droplet%position(:), droplet%velocity(:), velAir(:), droplet%radius)

        call search_refCELL(real(droplet%position(:)), droplet%refCellID)
        
    end subroutine

    subroutine output_mainDroplet(initial)
        logical, intent(in) :: initial
        character(255) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call mainDroplet%output_VTK(fname, deadline=initial)

        fname = case_dir//'/particle.csv'
        call mainDroplet%output_CSV(fname, dropletSolver%TimeOnSimu(step=timeStep, dimension=.true.), initial)

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

    subroutine divideAreaCoalescence_process
        ! integer, parameter :: num_divide = 4
        type(DropletGroup) dGroup
        integer i, j, k, id, m
        integer, allocatable :: ID_array(:)
        double precision AreaMin(3), AreaMax(3), width(3), delta(3), min_cdn(3), max_cdn(3)
        double precision, parameter :: deltaRatio = 1.d-2

        AreaMin(:) = 1.d9
        AreaMax(:) = -1.d9
        do m = 1, size(mainDroplet%droplet)
            if(mainDroplet%droplet(m)%status==0) then
                do i = 1, 3
                    AreaMin(i) = min(mainDroplet%droplet(m)%position(i), AreaMin(i))
                    AreaMax(i) = max(mainDroplet%droplet(m)%position(i), AreaMax(i))
                end do
            end if
        end do
        ! AreaMin(:) = AreaMin(:) - 1.d-9 ;print*, 'AreaMin:',AreaMin
        ! AreaMax(:) = AreaMax(:) + 1.d-9 ;print*, 'AreaMax:',AreaMax
        width(:) = (AreaMax(:) - AreaMin(:)) / dble(num_divide)   !;print*, 'width:',width
        delta(:) = (AreaMax(:) - AreaMin(:))*deltaRatio

        do k = 1, num_divide
            min_cdn(3) = AreaMin(3) + width(3)*dble(k-1) - delta(3)
            max_cdn(3) = AreaMin(3) + width(3)*dble(k) + delta(3)
            do j = 1, num_divide
                min_cdn(2) = AreaMin(2) + width(2)*dble(j-1) - delta(2)
                max_cdn(2) = AreaMin(2) + width(2)*dble(j) + delta(2)
                do i = 1, num_divide
                    min_cdn(1) = AreaMin(1) + width(1)*dble(i-1) - delta(1)
                    max_cdn(1) = AreaMin(1) + width(1)*dble(i) + delta(1)

                    ID_array = mainDroplet%IDinBox(min_cdn, max_cdn, status=0)
                    ! print*, 'divide_stat :', size(ID_array), min_cdn, max_cdn
                    dGroup%droplet = mainDroplet%droplet(ID_array)
                    call dGroup%coalescence_check()

                    do m = 1, size(ID_array)
                        id = ID_array(m)
                        mainDroplet%droplet(id) = dGroup%droplet(m)
                        if(dGroup%droplet(m)%coalesID > 0) mainDroplet%droplet(id)%coalesID = ID_array(dGroup%droplet(m)%coalesID)
                    end do

                end do
            end do
        end do

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
        print*, 'Now_Step_Time =', dropletSolver%TimeOnSimu(step=timeStep, dimension=.true.), '[sec]'
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

        TimeStart = dropletSolver%TimeOnSimu(step=n_start, dimension=.true.)
        TimeEnd = dropletSolver%TimeOnSimu(step=n_end, dimension=.true.)

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
            write(n_unit,'(A18, F18.2)') 'Temp [degC] =', dropletSolver%dropletEnvironment('Temperature')
            write(n_unit,'(A18, F18.2)') 'RH [%] =', dropletSolver%dropletEnvironment('RelativeHumidity')
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