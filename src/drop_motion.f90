module drop_motion_mod
    use flow_field
    use virusDroplet_m
    use equation_mod
    implicit none

    type, extends(virusDroplet_t) :: Droplet_onFlow
        integer :: adhesBoundID = 0
        type(reference_cell_t) refCELL
    end type Droplet_onFlow

    type(Droplet_onFlow), allocatable :: droplets(:)

    integer, private :: num_droplets   !全飛沫数
    integer interval
    real, private :: T, RH

    integer, target :: n_time  !時間ステップ
    integer, private :: num_restart
    integer n_start, n_end

    contains

    subroutine first_setting(case_dir)
        character(*), intent(in) :: case_dir

        call read_and_set_condition(case_dir)
        call set_coeff_drdt(T, RH)          !温湿度依存の係数の設定

        if(num_restart==0) then
            
            call random_set  !実行時刻に応じた乱数シード設定
            call calc_initial_position(case_dir)
            call calc_initial_radius
            call set_deathParam

        else if(num_restart==-1) then
            call read_initialDistribution(case_dir)

        else
            return  !リスタートなら無視

        end if

        call calc_minimumRadius(RH) !最小半径の計算

    end subroutine first_setting

    subroutine read_and_set_condition(dir)
        use filename_mod
        use path_operator_m
        character(*), intent(in) ::dir
        double precision dt, L, U
        double precision :: direction_g(3)
        character(99) path2FlowFile
        integer i, n_unit, num_drop

        OPEN(newunit=n_unit, FILE=dir//'/'//conditionFName, STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,*) num_restart
            read(n_unit,'()')
            read(n_unit,*) n_end
            read(n_unit,'()')
            read(n_unit,*) dt
            read(n_unit,'()')
            read(n_unit,*) interval
            read(n_unit,'()')
            read(n_unit,*) T
            read(n_unit,*) RH
            read(n_unit,'()')
            read(n_unit,*) num_drop
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

        if(num_restart == 0) call allocation_initialDroplets(num_drop) !通常実行なら割付

        if(num_restart > 0) then
            print*, 'Restart from', num_restart
        else if(num_restart == -1) then
            print*, 'InitialDistributon is Specified.'
        end if

        print*, 'n_end =',n_end
        print*, 'interval =',interval

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

    end subroutine read_and_set_condition

    subroutine set_initialDroplet(case_dir)
        character(*), intent(in) :: case_dir
        character(99) fname
        type(virusDroplet_t), allocatable :: droplets_ini(:)

        if(num_restart > 0) then

            print*, '**RESTART**'

            write(fname,'("'//case_dir//'/backup/backup_", i0, ".bu")') num_restart
            droplets = read_backup(fname)   !ここで自動割り付け

            n_start = num_restart
            n_time = n_start

        else
            droplets_ini = get_initialState_of_droplets()
            if(allocated(droplets)) deallocate(droplets)
            allocate(droplets(size(droplets_ini)))
            droplets(:)%virusDroplet_t = droplets_ini(:)
            n_start = 0
            n_time = n_start
            
        end if

        num_droplets = size(droplets)
        print*, 'num_droplets =', num_droplets
            
    end subroutine set_initialDroplet

    subroutine output_initialDroplet(case_dir)
        character(*), intent(in) :: case_dir

        if(num_restart <= 0) call output_droplet(case_dir, initial=.true.)  !リスタートでないなら初期配置出力

    end subroutine output_initialDroplet

    subroutine first_refCELLsearch
        integer j, num_drop
        logical success

        print*, 'first_refCELLsearch occured!'

        num_drop = size(droplets)

        if(unstructuredGrid) then
            j = 1
            droplets(j)%refCELL%ID = nearest_cell(real(droplets(j)%position(:)))

            droplets(j+1:)%refCELL%ID = droplets(j)%refCELL%ID !時間短縮を図る

            do j = 2, num_drop
                call search_refCELL(real(droplets(j)%position(:)), droplets(j)%refCELL%ID, stat=success)
                if(.not.success) droplets(j+1:)%refCELL%ID = droplets(j)%refCELL%ID
            end do

        else
            j = 1
            droplets(j)%refCELL%ID = get_cube_contains(real(droplets(j)%position(:)))    
            droplets(j)%refCELL%nodeID(:) = nearest_node(real(droplets(j)%position(:)), droplets(j)%refCELL%ID)

            droplets(j+1:)%refCELL = droplets(j)%refCELL !時間短縮を図る

            do j = 2, num_drop
                call search_refCELL_onCUBE(real(droplets(j)%position(:)), droplets(j)%refCELL)
            end do

        end if

    end subroutine first_refCELLsearch

    subroutine adhesion_check
        use adhesion_onSTL_m
        integer i
        
        if(unstructuredGrid) then
            do i = 1, num_droplets
                if(droplets(i)%status==0) then
                    call adhesion_onBound(droplets(i))
                    call area_check(droplets(i))
                end if
            end do
        else
            do i = 1, num_droplets
                if(droplets(i)%status==0) then
                    if(adhesion_onSTL(real(droplets(i)%position(:)))) call stop_droplet(droplets(i))
                    call area_check(droplets(i))
                end if
            end do
        end if

    end subroutine adhesion_check

    subroutine survival_check
        use terminalControler_m
        integer vfloat, vn
        ! double precision rand
        ! double precision, save :: death_rate = 0.d0
            
        vfloat = count(droplets(:)%status == 0)
        if(vfloat == 0) return  !浮遊数がゼロならリターン

        do vn = 1, num_droplets
            if ((droplets(vn)%status == 0).and.(droplets(vn)%deathParam > survival_rate(n_time))) then
                call stop_droplet(droplets(vn), status=-1)
            end if
        end do

        ! death_rate = death_rate + dble(vfloat)*(survival_rate(n_time) - survival_rate(n_time+1))    !このステップで死滅すべき飛沫数
        
        ! do while(death_rate >= 1.0d0)
        !     call random_number(rand)    !死滅IDは乱数で決まる
        !     vn = int(num_droplets*rand)
        !     if(vn < 1) cycle
        !     if (droplets(vn)%status == 0) then !浮遊粒子からのみ除去する
        !         droplets(vn)%status = -1
        !         ! droplets(vn)%position(:) = MIN_CDN(:) - 1.0d0 !計算エリア外に配置（不要かも）
        !         droplets(vn)%velocity(:) = 0.0d0
        !         death_rate = death_rate - 1.0d0

        !     end if
        ! end do
      
    end subroutine survival_check

    subroutine Calculation_Droplets
        implicit none
        integer vn

        if(unstructuredGrid) then
            !$omp parallel do
            do vn = 1, num_droplets
                if(droplets(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call evaporation(vn)    !蒸発方程式関連の処理
                call motionCalculation(vn)     !運動方程式関連の処理
            end do
            !$omp end parallel do

        else
            do vn = 1, num_droplets
                if(droplets(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call evaporation(vn)    !蒸発方程式関連の処理
                call motionCalculation_onCUBE(vn)     !運動方程式関連の処理
            end do

        end if

    end subroutine Calculation_Droplets

    subroutine evaporation(vn) !CALCULATE droplet evaporation
        integer, intent(in) :: vn
        double precision radius_n
      
        if (droplets(vn)%radius <= droplets(vn)%radius_min) return  !半径が最小になったものを除く
    
        radius_n = evaporatin_eq(droplets(vn)%radius)
        
        droplets(vn)%radius = max(radius_n, droplets(vn)%radius_min)
      
    end subroutine evaporation

    subroutine motionCalculation(vn)
        integer, intent(in) :: vn
        double precision velAir(3)

        velAir(:) = CELLs(droplets(vn)%refCELL%ID)%flowVelocity(:)
    
        call solve_motionEquation(droplets(vn)%position(:), droplets(vn)%velocity(:), velAir(:), droplets(vn)%radius)

        call search_refCELL(real(droplets(vn)%position(:)), droplets(vn)%refCELL%ID)
        
    end subroutine motionCalculation

    subroutine motionCalculation_onCUBE(vn)
        integer, intent(in) :: vn
        double precision velAir(3)
        type(reference_cell_t) :: RefC

        RefC = droplets(vn)%refCELL

        velAir(:) = get_velocity_f(RefC%nodeID, RefC%ID)

        call solve_motionEquation(droplets(vn)%position(:), droplets(vn)%velocity(:), velAir(:), droplets(vn)%radius)

        call search_refCELL_onCUBE(real(droplets(vn)%position(:)), droplets(vn)%refCELL)
    
    end subroutine motionCalculation_onCUBE
                    
    subroutine adhesion_onBound(droplet)
        use vector_m
        type(Droplet_onFlow), intent(inout) :: droplet
        integer JJ, JB, refCELL
        logical adhesion
        double precision :: r_vector(3), inner

        refCELL = droplet%refCELL%ID
        adhesion = .false.

        do JJ = 1, size(CELLs(refCELL)%boundFaceID)
            JB = CELLs(refCELL)%boundFaceID(JJ)

            r_vector(:) = droplet%position(:) - BoundFACEs(JB)%center(:)

            inner = dot_product(r_vector(:), BoundFACEs(JB)%normalVector(:))
            !外向き法線ベクトルと位置ベクトルの内積は、平面からの飛び出し量に相当

            if (inner + droplet%radius > 0.d0) then
                adhesion = .true.               !(飛び出し量+飛沫半径)がゼロ以上なら付着判定
                droplet%adhesBoundID = JB       !付着した境界面番号
            end if
        end do

        if (adhesion) call stop_droplet(droplet)

    end subroutine adhesion_onBound

    subroutine stop_droplet(droplet, status)
        type(Droplet_onFlow), intent(inout) :: droplet
        integer, optional :: status

        droplet%velocity(:) = 0.0d0
        if(present(status)) then
            droplet%status = status
        else
            droplet%status = 1
        end if

    end subroutine stop_droplet

    subroutine update_flow_check

        if(INTERVAL_FLOW <= 0) return

        call set_STEPinFLOW(Time_onSimulation(n_time))

        if(STEPinFLOW >= NextUpdate) call update_FlowField(first=.false.)   !流れ場の更新

    end subroutine update_flow_check

    subroutine update_FlowField(first)
        logical, intent(in) :: first

        if(INTERVAL_FLOW <= 0) then
            call read_steadyFlowData
        else
            if(first) call set_STEPinFLOW(Time_onSimulation(n_time))
            call read_unsteadyFlowData

        end if

        call set_MinMaxCDN

        if(first) then
            call preprocess_onFlowField         !流れ場の前処理
            if(n_start == 0) call first_refCELLsearch
        end if

        if(unstructuredGrid) then
            call boundary_setting(first)
            if(.not.first) call boundary_move
        end if
            
    end subroutine update_FlowField
                      
    subroutine boundary_move !境界面の移動に合わせて付着飛沫も移動
        integer vn, JB

        ! print*, 'CALL:boundary_move'

        do vn = 1, num_droplets
        
            if (droplets(vn)%status <= 0) cycle !付着していないならスルー
    
            JB = droplets(vn)%adhesBoundID
            if (JB > 0) then
                droplets(vn)%position(:) &
                    = droplets(vn)%position(:) + BoundFACEs(JB)%moveVector(:) !面重心の移動量と同じだけ移動
            else
                call area_check(droplets(vn))
    
            end if
        
        end do

        ! print*, 'FIN:boundary_move'

    end subroutine boundary_move

    subroutine area_check(droplet)
        type(Droplet_onFlow), intent(inout) :: droplet
        logical check
        integer L

        check = .false.
        do L = 1, 3
    
            if(droplet%position(L) < MIN_CDN(L)) then
                droplet%position(L) = MIN_CDN(L)
                check = .true.
            else if(droplet%position(L) > MAX_CDN(L)) then
                droplet%position(L) = MAX_CDN(L)
                check = .true.
            end if

        end do

        if (check) call stop_droplet(droplet)

    end subroutine area_check

    integer function drop_counter(name)
        character(*), intent(in) :: name

        select case(name)
            case('total')
                drop_counter = size(droplets)

            case('adhesion')
                drop_counter = count(droplets(:)%status > 0)

            case('floating')
                drop_counter = count(droplets(:)%status == 0)

            case('death')
                drop_counter = count(droplets(:)%status == -1)

            case('coalescence')
                drop_counter = count(droplets(:)%status == -2)

            case default
                print*, '**ERROR [drop_counter] : ', name, ' is not found.**'
                stop

        end select

    end function drop_counter

    subroutine coalescence_check
        use terminalControler_m
        integer d1, d2, floatings
        integer, save :: last_coalescence = 0, last_floatings = 0
        double precision :: distance, r1, r2

        floatings = drop_counter('floating')
        if(floatings > last_floatings) last_coalescence = n_time    !浮遊数が増加したら付着判定再起動のため更新
        last_floatings = floatings

        !最後の合体から100ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((n_time - last_coalescence) > 100) return

        call set_formatTC('(" Coalescence_check [step:", i10, "/", i10, "]")')
        call print_sameLine([n_time, last_coalescence+100])

        !$OMP parallel do private(distance, r1, r2)
        drop1 : do d1 = 1, num_droplets - 1
            if(droplets(d1)%status/=0) cycle drop1

            r1 = droplets(d1)%radius

            drop2 : do d2 = d1 + 1, num_droplets
                if(droplets(d2)%status/=0) cycle drop2

                r2 = droplets(d2)%radius

                distance = norm2(droplets(d2)%position(:) - droplets(d1)%position(:))

                if((r1+r2) >= distance) then
                    ! print*, d1, 'and', d2, 'coalesce!'
                    if(r1 >= r2) then
                        call coalescence(droplets(d1), droplets(d2))
                    else
                        call coalescence(droplets(d2), droplets(d1))
                    end if
                    last_coalescence = n_time

                end if

            end do drop2

        end do drop1
        !$OMP end parallel do

    end subroutine coalescence_check

    subroutine coalescence(droplet1, droplet2)
        type(Droplet_onFlow), intent(inout) :: droplet1, droplet2
        double precision :: radius_c, volume1, volume2, velocity_c(3)

        volume1 = droplet1%radius**3
        volume2 = droplet2%radius**3
        radius_c = (volume1 + volume2)**(1.d0/3.d0) !結合後の飛沫半径
        velocity_c(:) = (volume1*droplet1%velocity(:) + volume2*droplet2%velocity(:)) / (volume1 + volume2)
        
        droplet1%radius = radius_c
        droplet1%velocity(:) = velocity_c(:)

        droplet2%radius = 0.d0
        call stop_droplet(droplet2, status=-2)
        
    end subroutine coalescence

    function read_backup(fname) result(droplets_read)
        implicit none
        character(*), intent(in) :: fname
        type(Droplet_onFlow), allocatable :: droplets_read(:)
        integer i, n_unit, num_drop
    
        print*, 'READ:', trim(fname)
        open(newunit=n_unit, form='unformatted', file=fname, status='old')
            read(n_unit) num_drop

            allocate(droplets_read(num_drop))

            do i = 1, num_drop
                read(n_unit) droplets_read(i)
            end do
        close(n_unit)
    
      
    end function read_backup

    subroutine output_backup(dir)
        character(*), intent(in) :: dir
        integer i, n_unit, num_drop
        character(99) fname

        num_drop = size(droplets(:))
        write(fname,'("'//dir//'/backup_", i0, ".bu")') n_time

        open(newunit=n_unit, form='unformatted', file=fname, status='replace')
            write(n_unit) num_drop
            do i = 1, num_drop
                write(n_unit) droplets(i)
            end do
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine output_backup

    subroutine output_droplet(case_dir, initial)
        character(*), intent(in) :: case_dir
        logical, intent(in) :: initial
        character(99) fname

        write(fname,'("'//case_dir//'/VTK/drop_", i0, ".vtk")') n_time
        call output_droplet_VTK(fname, droplets(:)%virusDroplet_t, initial)

        fname = case_dir//'/particle.csv'
        call output_droplet_CSV(fname, droplets(:)%virusDroplet_t, Time_onSimulation(n_time, dimension=.true.), initial)

        if(.not.initial) call output_backup(case_dir//'/backup')

    end subroutine

    real function environment(name)
        character(*), intent(in) :: name

        select case(name)
            case('Temperature')
                environment = T
            case('RelativeHumidity')
                environment = RH
            case default
                environment = -1.e20
        end select

    end function environment

end module drop_motion_mod