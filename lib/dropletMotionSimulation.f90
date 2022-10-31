module dropletMotionSimulation
    use virusDroplet_m
    use dropletGenerator_m
    use dropletEquation_m
    use flow_field_m
    use timeKeeper_m
    implicit none

    private

    integer outputInterval

    integer num_restart
    integer, target :: timeStep
    integer n_start, n_end

    type(TimeKeeper) tK

    character(:), allocatable :: case_dir

    logical :: startFlag = .false.
    integer :: last_coalescenceStep=0
    logical generationFlag

    logical :: adhesionSwitch = .true.
    integer :: coalescenceLimit=10000, num_divide=4
    character(:), allocatable :: radiusDistributionFilename

    type(virusDroplet_t), allocatable :: mainDroplets(:)

    type(DropletEquationSolver), target :: dropletSolver

    type(DropletGenerator) dropGenerator

    type(FlowField) flow_field

    public mainDropletLoop, simulationSetUp, output_ResultSummary, read_basicSettingOnSimulation

    contains

    subroutine simulationSetUp(case_name)
        use virusDroplet_m
        use conditionValue_m
        character(*), intent(in) :: case_name
        type(conditionValue_t) condVal

        case_dir = case_name
        call create_CaseDirectory

        condVal = read_condition(case_dir)

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

        timeStep = n_start

        dropGenerator = DropletGenerator_( &
                            dropletSolver, radiusDistributionFilename, case_dir, &
                            generationRate = condVal%periodicGeneration &
                        )

        if(num_restart <= 0) then

            if(condVal%isInitialDistributionSpecified()) then

                mainDroplets = read_backup(case_dir//'/'//condVal%initialDistributionFName)

            else

                mainDroplets = dropGenerator%generateDroplet(condVal%num_drop, TimeOnSimu())

            end if

            call output_mainDroplet(initial = .true.)   !この時点では、飛沫の参照セルは見つかっていない

        else

            print*, '**RESTART**'

            block
                character(255) fname

                write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') num_restart
                mainDroplets = read_backup(trim(fname))   !ここで自動割り付け

            end block

        end if

        print*, 'num_droplets =', size(mainDroplets)
        last_coalescenceStep = 0

        call checkpoint

        if(condVal%isMeshFileSpecified()) then
            flow_field = FlowField_( &                !流れ場の取得
                dropletSolver%TimeStep2RealTime(step=n_start, dimension=.false.), &
                condVal%PATH2FlowFile, condVal%DT_FLOW, condVal%OFFSET, condVal%INTERVAL_FLOW, &
                condVal%LoopHead, condVal%LoopTail, &
                condVal%meshFile &
            )

        else
            flow_field = FlowField_( &                !流れ場の取得
                dropletSolver%TimeStep2RealTime(step=n_start, dimension=.false.), &
                condVal%PATH2FlowFile, condVal%DT_FLOW, condVal%OFFSET, condVal%INTERVAL_FLOW, &
                condVal%LoopHead, condVal%LoopTail &
            )

        end if
            
        if(num_restart <= 0) call first_refCellSearch(mainDroplets)

    end subroutine

    subroutine read_basicSettingOnSimulation
        integer n_unit
        character(23) :: fname ='option/basicSetting.nml'
        character(255) radiusDistributionFNAME
        namelist /basicSetting/ coalescenceLimit, adhesionSwitch, num_divide, radiusDistributionFNAME

        print*, 'READ : ', fname
        
        open(newunit=n_unit, file=fname, status='old', action='read')
            read(n_unit, nml=basicSetting)
        close(n_unit)

        radiusDistributionFilename = trim(radiusDistributionFNAME)

    end subroutine

    subroutine mainDropletLoop
        integer, pointer :: n => timeStep
        
        print '("*******************************************")'
        print '("            START step_loop                ")'
        print '("*******************************************")'

        do n = n_start + 1, n_end           !ステップ数だけループ
            
            call dropGenerator%periodicGeneration(mainDroplets, TimeOnSimu(), generationFlag)

            if(adhesionSwitch) call adhesion_check(mainDroplets)

            call survival_check(mainDroplets, TimeOnSimu())           !生存率に関する処理

            call coalescence_process        !飛沫間の合体判定

            call Calculation_Droplets     !飛沫の運動計算

            if (mod(n, outputInterval) == 0) call periodicOutput(n)            !出力

            call check_FlowFieldUpdate        !流れ場の更新チェック

        end do

        print '("*******************************************")'
        print '("             END step_loop                 ")'
        print '("*******************************************")'

    end subroutine

    subroutine check_FlowFieldUpdate

        if(timeStep==n_end) return

        call flow_field%set_time(TimeOnSimu())

        if(flow_field%isUpdateTiming()) then
            call flow_field%update()   !流れ場の更新
            call dropletOnBoundary(mainDroplets)
        end if

    end subroutine

    subroutine first_refCellSearch(droplets)
        type(virusDroplet_t), intent(inout) :: droplets(:)
        integer j, num_drop
        logical success
        real position(3)

        print*, 'first_refCellSearch occured!'

        num_drop = size(droplets)

        j = 1
        droplets(j)%refCellID = flow_field%nearest_cell(real(droplets(j)%position(:)))

        droplets(j+1:)%refCellID = droplets(j)%refCellID !時間短縮を図る

        do j = 2, num_drop
            position = real(droplets(j)%position(:))
            call flow_field%search_refCELL(position, droplets(j)%refCellID, stat=success)
            if(.not.success) droplets(j+1:)%refCellID = droplets(j)%refCellID
        end do
        
        ! call mainMesh%sort()

    end subroutine
   
    subroutine adhesion_check(droplets)
        use unstructuredGrid_m
        type(virusDroplet_t), intent(inout) :: droplets(:)
        integer i

        do i = 1, size(droplets)
            if(droplets(i)%isFloating()) then
                call flow_field%adhesionCheckOnBound( &
                    droplets(i)%position, droplets(i)%get_radius(), droplets(i)%refCellID, &
                    stat=droplets(i)%adhesBoundID &
                    )
                if (droplets(i)%adhesBoundID >= 1) call droplets(i)%stop_droplet()
            end if
        end do
            
        call area_check(droplets)

    end subroutine

    subroutine area_check(droplets)
        type(virusDroplet_t), intent(inout) :: droplets(:)
        logical check
        real areaMin(3), areaMax(3)
        integer i, J

        call flow_field%get_MinMaxOfGrid(areaMin, areaMax)

        do i = 1, size(droplets)
            check = .false.
            do J = 1, 3
        
                if(droplets(i)%position(J) < areaMin(J)) then
                    droplets(i)%position(J) = areaMin(J)
                    check = .true.
                else if(droplets(i)%position(J) > areaMax(J)) then
                    droplets(i)%position(J) = areaMax(J)
                    check = .true.
                end if

            end do

            if (check) call droplets(i)%stop_droplet()

        end do

    end subroutine
                      
    subroutine dropletOnBoundary(droplets) !境界面の移動に合わせて付着飛沫も移動
        use unstructuredGrid_m
        type(virusDroplet_t), intent(inout) :: droplets(:)
        integer vn, JB

        ! print*, 'CALL:dropletOnBoundary'

        do vn = 1, size(droplets)
        
            if (droplets(vn)%isFloating()) cycle !付着していないならスルー
    
            JB = droplets(vn)%adhesBoundID
            if (JB > 0) then
                droplets(vn)%position(:) &
                    = droplets(vn)%position(:) &
                    + flow_field%get_movementVectorOfBoundarySurface(JB) !面重心の移動量と同じだけ移動
            end if
        
        end do

        ! call area_check(dGroup)

        ! print*, 'FIN:dropletOnBoundary'

    end subroutine

    subroutine Calculation_Droplets()
        integer vn, targetID

        !$omp parallel do
        do vn = 1, size(mainDroplets)

            if(mainDroplets(vn)%isFloating())then

                call evaporationProcess(mainDroplets(vn))    !蒸発方程式関連の処理
                call motionCalculation(mainDroplets(vn))     !運動方程式関連の処理
            
            else 
                targetID = mainDroplets(vn)%coalescenceID()
                if(targetID > 0) then
                    !合体飛沫の片割れも移動させる
                    mainDroplets(vn)%position = mainDroplets(targetID)%position
                    mainDroplets(vn)%velocity = mainDroplets(targetID)%velocity
                end if
            end if

        end do
        !$omp end parallel do

    end subroutine

    subroutine evaporationProcess(droplet) !CALCULATE droplet evaporation
        use virusDroplet_m
        type(virusDroplet_t) droplet
        double precision dr
      
        if (.not.droplet%isEvaporating()) return  !半径が最小になったものを除く
    
        dr = dropletSolver%evaporationEq(droplet%get_radius())  !半径変化量
        
        call droplet%evaporation(dr)
      
    end subroutine


    subroutine motionCalculation(droplet)
        use virusDroplet_m
        type(virusDroplet_t) droplet
        double precision velAir(3)
        real position(3)

        velAir(:) = flow_field%get_flowVelocityInCELL(droplet%refCellID)
    
        call dropletSolver%solve_motionEquation(droplet%position(:), droplet%velocity(:), velAir(:), droplet%get_radius())

        position = real(droplet%position(:))
        call flow_field%search_refCELL(position, droplet%refCellID)
        
    end subroutine

    subroutine output_mainDroplet(initial)
        use filename_m, only : IniDistributionFName => InitialDistributionFileName
        logical, intent(in) :: initial
        character(255) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call output_droplet_VTK(mainDroplets, fname, deadline=initial)

        fname = case_dir//'/particle.csv'
        call output_droplet_CSV(mainDroplets, fname, TimeOnSimu(dimension=.true.), initial)

        if(initial) then
            fname = case_dir//'/backup/'//IniDistributionFName
        else
            write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') timeStep
        end if
        call output_backup(mainDroplets, trim(fname))

    end subroutine

    subroutine coalescence_process
        use terminalControler_m
        integer numFloating, num_coalescence
        
        numFloating = dropletCounter(mainDroplets, 'floating')
        if(generationFlag) last_coalescenceStep = timeStep - 1    !飛沫発生が起こったら前ステップに付着が起こったことにする（付着判定再起動のため）

        !最後の合体から指定ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((timeStep - last_coalescenceStep) > coalescenceLimit) return

        call set_formatTC('(" Coalescence_check [step:", i10, "/", i10, "]")')
        call print_progress([timeStep, last_coalescenceStep + coalescenceLimit])

        ! call mainDroplet%coalescence_check(stat = num_coalescence)
        call divideAreaCoalescence_process(num_coales = num_coalescence)

        if(num_coalescence >= 1) last_coalescenceStep = timeStep

    end subroutine

    subroutine divideAreaCoalescence_process(num_coales)
        type(virusDroplet_t), allocatable :: droplets(:)
        integer, intent(out) :: num_coales
        integer i, j, k, id, m, stat_coales
        integer, allocatable :: ID_array(:)
        double precision AreaMin(3), AreaMax(3), width(3), delta(3), min_cdn(3), max_cdn(3)
        double precision, parameter :: deltaRatio = 1.d-2

        num_coales = 0
        if(num_divide <= 0) return

        call get_dropletsArea(mainDroplets, AreaMin, AreaMax)

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

                    ID_array = dropletIDinBox(mainDroplets, min_cdn, max_cdn, status=0)
                    ! print*, 'divide_stat :', size(ID_array), min_cdn, max_cdn
                    droplets = mainDroplets(ID_array)  !分割エリア内の飛沫を抽出（ここでIDが変わる）
                    call coalescence_check(droplets, stat=stat_coales)  !分割エリア内で合体判定
                    num_coales = num_coales + stat_coales

                    block
                        integer coalesID
                        do m = 1, size(ID_array)
                            id = ID_array(m)
                            mainDroplets(id) = droplets(m)  !飛沫情報をもとのIDに格納

                            coalesID = droplets(m)%coalescenceID()
                            !合体飛沫については、合体先ID（coalesID）ももとのIDに戻す必要がある
                            if(coalesID > 0) mainDroplets(id)%coalesID = ID_array(coalesID)
                        end do
                    end block

                end do
            end do
        end do

    end subroutine

    subroutine checkpoint
        character(1) input

        do while(.not.startFlag)
            print*, 'Do you want to start the calculation? (y/n)'
            read(5,*) input

            select case(input)
            case('y')
                startFlag = .true.

            case('n')
                stop

            end select

        end do

        tk = TimeKeeper_()

    end subroutine

    subroutine create_CaseDirectory
        use path_operator_m

        call make_directory(case_dir//'/VTK')
        call make_directory(case_dir//'/backup')
        
    end subroutine

    subroutine periodicOutput(nowStep)
        use terminalControler_m
        integer, intent(in) :: nowStep
        real nearerSearchFalseRate

        print*, '[Start Date] ' // tk%startDateAndTime()
        print '(" ** It will take", f8.2, " minites **")', real(n_end - nowStep)/(60.*real(nowStep)/tK%erapsedTime())
        print*, 'Now_Step_Time =', TimeOnSimu(dimension=.true.), '[sec]'
        print*, '# floating :', dropletCounter(mainDroplets, 'floating'), '/', dropletCounter(mainDroplets, 'total')
        nearerSearchFalseRate = flow_field%get_nearerSearchFalseRate()
        if(nearerSearchFalseRate >= 1.) print*, '# searchFalseRate :', nearerSearchFalseRate, '%'
        call output_mainDroplet(initial=.false.)
        print '("====================================================")'
        call reset_formatTC

    end subroutine

    subroutine output_ResultSummary()
        use dropletEquation_m
        integer n_unit, cnt
        real erapsed_time
        character(50) fname
        logical existance
        double precision TimeStart, TimeEnd
        character(:), allocatable :: startDAT, endDAT
        
        erapsed_time = tk%erapsedTime()

        startDAT = '[Start Date] ' // tk%startDateAndTime()
        endDAT = '[ END  Date] ' // nowDateAndTime()
        print*, startDAT
        print*, endDAT

        fname = case_dir//'/ResultSummary.txt'
        inquire(file=fname, exist=existance)
        cnt = 0
        do while(existance)
            cnt = cnt + 1
            write(fname,'("'//case_dir//'/ResultSummary_", i0, ".txt")') cnt
            inquire(file=fname, exist=existance)
        end do

        TimeStart = dropletSolver%TimeStep2RealTime(step=n_start, dimension=.true.)
        TimeEnd = dropletSolver%TimeStep2RealTime(step=n_end, dimension=.true.)

        open(newunit=n_unit, file=fname, status='new')
            write(n_unit,*)'*******************************************'
            write(n_unit,*)'*                                         *'
            write(n_unit,*)'*             Result Summary              *'
            write(n_unit,*)'*                                         *'
            write(n_unit,*)'*******************************************'
            write(n_unit,'(A)') '======================================================='
            write(n_unit,*) startDAT
            write(n_unit,*) endDAT
            write(n_unit, '(A18, F15.3, 2X, A)') 'Erapsed Time =', erapsed_time, '[sec]'
            write(n_unit, '(A18, F15.3, 2X, A)') 'Cost of Calc =', &
                    erapsed_time / (TimeEnd - TimeStart), '[sec/sec]'
            write(n_unit,'(A)') '======================================================='
            write(n_unit, '(A18, 2(F15.3,2x,A))') 'Time [sec] =', TimeStart, '-', TimeEnd
            write(n_unit, '(A18, 2(I15,2x,A))') 'Step =', n_start, '-', n_end !計算回数
            write(n_unit,'(A18, I15)') 'OutputInterval =', outputInterval
            write(n_unit,'(A)') '======================================================='
            write(n_unit,'(A18, I15)') '#Droplets =', dropletCounter(mainDroplets, 'total')
            write(n_unit,'(A18, I15)') 'floating =', dropletCounter(mainDroplets, 'floating')
            write(n_unit,'(A18, I15)') 'death =', dropletCounter(mainDroplets, 'death') !生存率で消滅
            write(n_unit,'(A18, I15)') 'coalescence =', dropletCounter(mainDroplets, 'coalescence') !生存率で消滅
            write(n_unit,'(A18, I15)') 'adhesion =', dropletCounter(mainDroplets, 'adhesion') !付着したすべてのウイルス数
            write(n_unit,'(A)') '======================================================='
            write(n_unit,'(A18, F18.2)') 'Temp [degC] =', dropletSolver%dropletEnvironment('Temperature')
            write(n_unit,'(A18, F18.2)') 'RH [%] =', dropletSolver%dropletEnvironment('RelativeHumidity')
            write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', flow_field%get_defaultFlowFileName()
            write(n_unit, '(A18, I15,2x,A, f8.3,2x,A)') 'SearchFalseInfo :', flow_field%get_num_nearerSearchFalse(), &
                    ' (', flow_field%get_nearerSearchFalseRate(), '%)'

        close(n_unit)
        
    end subroutine

    double precision function TimeOnSimu(dimension)
        logical, intent(in), optional :: dimension

        if(present(dimension)) then
            TimeOnSimu = dropletSolver%TimeStep2RealTime(timeStep, dimension)
        else
            TimeOnSimu = dropletSolver%TimeStep2RealTime(timeStep, .false.)
        end if

    end function

    function DateAndTime_string(date, time) result(string)
        character(*), intent(in) :: date, time
        character(:), allocatable :: string

        string = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' ' &
              //time(1:2)//':'//time(3:4)//':'//time(5:6)

    end function

    ! subroutine random_set    !実行時刻に依存した乱数シードを指定する
    !     implicit none
    !     integer :: seedsize, i
    !     integer, allocatable :: seed(:)

    !     print*, 'call:random_set'
    
    !     call random_seed(size=seedsize) !シードのサイズを取得。（コンパイラごとに異なるらしい）
    !     allocate(seed(seedsize)) !新シード配列サイズの割り当て
    
    !     do i = 1, seedsize
    !         call system_clock(count=seed(i)) !時間を新シード配列に取得
    !     end do
    
    !     call random_seed(put=seed(:)) !新シードを指定
          
    ! end subroutine random_set

end module dropletMotionSimulation