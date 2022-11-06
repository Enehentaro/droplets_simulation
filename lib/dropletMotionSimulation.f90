module dropletMotionSimulation
    use virusDroplet_m
    use dropletGenerator_m
    use dropletEquation_m
    use flow_field_m
    use timeKeeper_m
    implicit none

    private

    integer outputInterval

    type(TimeKeeper) tK

    character(:), allocatable :: case_dir

    logical :: startFlag = .false.
    integer :: last_coalescenceStep=0
    logical generationFlag

    logical :: adhesionSwitch = .true.
    integer :: coalescenceLimit=10000, num_divide=4
    character(:), allocatable :: radiusDistributionFilename

    public RunDropletsSimulation, read_basicSettingOnSimulation

    contains

    subroutine simulationSetUp(case_name, droplets, dropletSolver, dropGenerator, flow_field, n_start, n_end)
        use virusDroplet_m
        use conditionValue_m
        character(*), intent(in) :: case_name
        type(virusDroplet_t), allocatable, intent(out) :: droplets(:)
        type(DropletEquationSolver), target, intent(out) :: dropletSolver
        type(DropletGenerator), intent(out) :: dropGenerator
        type(FlowField), intent(out) :: flow_field
        integer, intent(out) :: n_start, n_end
        type(conditionValue_t) condVal

        case_dir = case_name
        call create_CaseDirectory

        condVal = read_condition(case_dir)

        n_end = condVal%stepEnd
        outputInterval = condVal%outputInterval
        print*, 'n_end =',n_end
        print*, 'output interval =', outputInterval

        dropletSolver = DropletEquationSolver_( &
                            condVal%dt, condVal%L, condVal%U, &
                            condVal%direction_g, condVal%T, condVal%RH &
                        )

        n_start = max(condVal%restart, 0)

        dropGenerator = DropletGenerator_( &
                            dropletSolver, radiusDistributionFilename, case_dir, &
                            generationRate = condVal%periodicGeneration &
                        )

        if(condVal%restart <= 0) then

            if(condVal%isInitialDistributionSpecified()) then

                droplets = read_backup(case_dir//'/'//condVal%initialDistributionFName)

            else

                droplets = dropGenerator%generateDroplet(condVal%num_drop, dropletSolver%TimeStep2RealTime(n_start, .false.))

            end if

            call output_droplet_process(initial=.true., droplets=droplets, &
                timeStep=n_start, real_time=dropletSolver%TimeStep2RealTime(n_start, .true.))   !この時点では、飛沫の参照セルは見つかっていない

        else

            print*, '**RESTART**'

            block
                character(255) fname

                write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') condVal%restart
                droplets = read_backup(trim(fname))   !ここで自動割り付け

            end block

        end if

        print*, 'num_droplets =', size(droplets)
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
            
        if(condVal%restart <= 0) call first_refCellSearch(droplets, flow_field)

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

    subroutine RunDropletsSimulation(case_name)
        character(*), intent(in) :: case_name
        type(virusDroplet_t), allocatable :: mainDroplets(:)
        type(DropletEquationSolver), target :: dropletSolver
        type(DropletGenerator) dropGenerator
        type(FlowField) flow_field
        integer n, n_start, n_end
        double precision timeInSimulation

        call simulationSetUp(case_name, mainDroplets, dropletSolver, dropGenerator, flow_field, n_start, n_end)
        
        print '("*******************************************")'
        print '("            START step_loop                ")'
        print '("*******************************************")'

        do n = n_start + 1, n_end           !ステップ数だけループ
            timeInSimulation = dropletSolver%TimeStep2RealTime(n, .false.)
            
            call dropGenerator%periodicGeneration(mainDroplets, timeInSimulation, generationFlag)

            if(adhesionSwitch) call adhesion_check(mainDroplets, flow_field)

            call survival_check(mainDroplets, timeInSimulation)           !生存率に関する処理

            call coalescence_process(mainDroplets, n)        !飛沫間の合体判定

            call Calculation_Droplets(mainDroplets, dropletSolver, flow_field)     !飛沫の運動計算

            if (mod(n, outputInterval) == 0) call periodicOutput(n, n_end, mainDroplets, &
                dropletSolver%TimeStep2RealTime(n, .true.),  flow_field%get_nearerSearchFalseRate())!出力

            if(n==n_end) exit
            call check_FlowFieldUpdate(flow_field, mainDroplets, timeInSimulation)        !流れ場の更新チェック

        end do

        print '("*******************************************")'
        print '("             END step_loop                 ")'
        print '("*******************************************")'

        call output_ResultSummary(mainDroplets, dropletSolver, flow_field, n_start, n_end)

    end subroutine

    subroutine check_FlowFieldUpdate(flow_field, droplets, time)
        type(FlowField), intent(inout) ::  flow_field
        type(virusDroplet_t), intent(inout) ::  droplets(:)
        double precision, intent(in) :: time

        call flow_field%set_time(time)

        if(flow_field%isUpdateTiming()) then
            call flow_field%update()   !流れ場の更新
            call dropletOnBoundary(droplets, flow_field)
        end if

    end subroutine

    subroutine first_refCellSearch(droplets, flow_field)
        type(virusDroplet_t), intent(inout) :: droplets(:)
        type(FlowField), intent(inout) ::  flow_field
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
   
    subroutine adhesion_check(droplets, flow_field)
        use unstructuredGrid_m
        type(virusDroplet_t), intent(inout) :: droplets(:)
        type(FlowField), intent(in) ::  flow_field
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
            
        call area_check(droplets, flow_field)

    end subroutine

    subroutine area_check(droplets, flow_field)
        type(virusDroplet_t), intent(inout) :: droplets(:)
        type(FlowField), intent(in) ::  flow_field
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
                      
    subroutine dropletOnBoundary(droplets, flow_field) !境界面の移動に合わせて付着飛沫も移動
        use unstructuredGrid_m
        type(virusDroplet_t), intent(inout) :: droplets(:)
        type(FlowField), intent(in) ::  flow_field
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

    subroutine Calculation_Droplets(droplets, dropletSolver, flow_field)
        type(virusDroplet_t), intent(inout) :: droplets(:)
        type(DropletEquationSolver), intent(in) :: dropletSolver
        type(FlowField), intent(in) ::  flow_field
        integer vn, targetID
        double precision velAir(3)
        real position(3)

        !$omp parallel do
        do vn = 1, size(droplets)

            if(droplets(vn)%isFloating())then

                call evaporationProcess(droplets(vn), dropletSolver)    !蒸発方程式関連の処理
               
                velAir(:) = flow_field%get_flowVelocityInCELL(droplets(vn)%refCellID)
                call dropletSolver%solve_motionEquation(&
                    droplets(vn)%position(:), droplets(vn)%velocity(:), velAir(:), droplets(vn)%get_radius())
        
                position = real(droplets(vn)%position(:))
                call flow_field%search_refCELL(position, droplets(vn)%refCellID)
            
            else 
                targetID = droplets(vn)%coalescenceID()
                if(targetID > 0) then
                    !合体飛沫の片割れも移動させる
                    droplets(vn)%position = droplets(targetID)%position
                    droplets(vn)%velocity = droplets(targetID)%velocity
                end if
            end if

        end do
        !$omp end parallel do

    end subroutine

    subroutine evaporationProcess(droplet, dropletSolver) !CALCULATE droplet evaporation
        use virusDroplet_m
        type(virusDroplet_t), intent(inout) :: droplet
        type(DropletEquationSolver), intent(in) :: dropletSolver
        double precision dr
      
        if (.not.droplet%isEvaporating()) return  !半径が最小になったものを除く
    
        dr = dropletSolver%evaporationEq(droplet%get_radius())  !半径変化量
        
        call droplet%evaporation(dr)
      
    end subroutine

    subroutine output_droplet_process(initial, droplets, timeStep, real_time)
        use filename_m, only : IniDistributionFName => InitialDistributionFileName
        logical, intent(in) :: initial
        type(virusDroplet_t), intent(in) :: droplets(:)
        integer, intent(in) :: timeStep
        double precision, intent(in) :: real_time

        character(255) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') timeStep
        call output_droplet_VTK(droplets, fname, deadline=initial)

        fname = case_dir//'/particle.csv'
        call output_droplet_CSV(droplets, fname, real_time, initial)

        if(initial) then
            fname = case_dir//'/backup/'//IniDistributionFName
        else
            write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') timeStep
        end if
        call output_backup(droplets, trim(fname))

    end subroutine

    subroutine coalescence_process(mainDroplets, timeStep)
        use terminalControler_m
        type(virusDroplet_t), intent(inout) :: mainDroplets(:)
        integer, intent(in) :: timeStep
        integer numFloating, num_coalescence
        
        numFloating = dropletCounter(mainDroplets, 'floating')
        if(generationFlag) last_coalescenceStep = timeStep - 1    !飛沫発生が起こったら前ステップに付着が起こったことにする（付着判定再起動のため）

        !最後の合体から指定ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((timeStep - last_coalescenceStep) > coalescenceLimit) return

        call set_formatTC('(" Coalescence_check [step:", i10, "/", i10, "]")')
        call print_progress([timeStep, last_coalescenceStep + coalescenceLimit])

        ! call mainDroplet%coalescence_check(stat = num_coalescence)
        call divideAreaCoalescence_process(num_coales = num_coalescence, mainDroplets = mainDroplets)

        if(num_coalescence >= 1) last_coalescenceStep = timeStep

    end subroutine

    subroutine divideAreaCoalescence_process(num_coales, mainDroplets)
        integer, intent(out) :: num_coales
        type(virusDroplet_t), intent(inout) :: mainDroplets(:)
        type(virusDroplet_t), allocatable :: droplets(:)
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

    subroutine periodicOutput(nowStep, endStep, droplets, real_time, nearerSearchFalseRate)
        use terminalControler_m
        integer, intent(in) :: nowStep, endStep
        type(virusDroplet_t), intent(in) :: droplets(:)
        double precision, intent(in) :: real_time
        real, intent(in) ::  nearerSearchFalseRate

        print*, '[Start Date] ' // tk%startDateAndTime()
        print '(" ** It will take", f8.2, " minites **")', real(endStep - nowStep)/(60.*real(nowStep)/tK%erapsedTime())
        print*, 'Now_Step_Time =', real_time, '[sec]'
        print*, '# floating :', dropletCounter(droplets, 'floating'), '/', dropletCounter(droplets, 'total')
        if(nearerSearchFalseRate >= 1.) print*, '# searchFalseRate :', nearerSearchFalseRate, '%'
        call output_droplet_process(initial=.false., droplets=droplets, timeStep=nowStep, real_time=real_time)
        print '("====================================================")'
        call reset_formatTC

    end subroutine

    subroutine output_ResultSummary(droplets, dropletSolver, flow_field, n_start, n_end)
        use dropletEquation_m
        type(virusDroplet_t), allocatable, intent(in) :: droplets(:)
        type(DropletEquationSolver), target, intent(in) :: dropletSolver
        type(FlowField), intent(in) :: flow_field
        integer, intent(in) :: n_start, n_end
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
            write(n_unit,'(A18, I15)') '#Droplets =', dropletCounter(droplets, 'total')
            write(n_unit,'(A18, I15)') 'floating =', dropletCounter(droplets, 'floating')
            write(n_unit,'(A18, I15)') 'death =', dropletCounter(droplets, 'death') !生存率で消滅
            write(n_unit,'(A18, I15)') 'coalescence =', dropletCounter(droplets, 'coalescence') !生存率で消滅
            write(n_unit,'(A18, I15)') 'adhesion =', dropletCounter(droplets, 'adhesion') !付着したすべてのウイルス数
            write(n_unit,'(A)') '======================================================='
            write(n_unit,'(A18, F18.2)') 'Temp [degC] =', dropletSolver%dropletEnvironment('Temperature')
            write(n_unit,'(A18, F18.2)') 'RH [%] =', dropletSolver%dropletEnvironment('RelativeHumidity')
            write(n_unit,'(A18, 2X, A)') 'Used FlowFile :', flow_field%get_defaultFlowFileName()
            write(n_unit, '(A18, I15,2x,A, f8.3,2x,A)') 'SearchFalseInfo :', flow_field%get_num_nearerSearchFalse(), &
                    ' (', flow_field%get_nearerSearchFalseRate(), '%)'

        close(n_unit)
        
    end subroutine

    ! double precision function TimeOnSimu(dropletSolver, dimension)
    !     type(DropletEquationSolver), target, intent(in) :: dropletSolver
    !     logical, intent(in), optional :: dimension

    !     if(present(dimension)) then
    !         TimeOnSimu = dropletSolver%TimeStep2RealTime(timeStep, dimension)
    !     else
    !         TimeOnSimu = dropletSolver%TimeStep2RealTime(timeStep, .false.)
    !     end if

    ! end function

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