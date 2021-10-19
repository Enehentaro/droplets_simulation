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

    type pathCase_t
        character(:), allocatable :: DIR, VTK, backup
    end type pathCase_t

    type(pathCase_t) path

    integer, target :: n_time  !時間ステップ
    integer, private :: num_restart
    integer n_start, n_end
    double precision, private :: Rdt    !飛沫計算と気流計算の時間間隔の比

    contains

    subroutine first_setting

        call read_and_set_condition
        call set_coeff_drdt(T, RH)          !温湿度依存の係数の設定

        if(num_restart==0) then
            
            call random_set  !実行時刻に応じた乱数シード設定
            call calc_initial_position(path%DIR)
            call calc_initial_radius
            call set_death_param

        else if(num_restart==-1) then
            call read_initialDistribution

        else
            return  !リスタートなら無視

        end if

        call calc_minimumRadius(RH) !最小半径の計算

    end subroutine first_setting

    subroutine set_case_path(case_name)
        character(*), intent(in) :: case_name
        path%DIR = trim(case_name)//'/'
        path%VTK = trim(case_name)//'/VTK/'
        path%backup = trim(case_name)//'/backup/'
    end subroutine set_case_path

    subroutine read_and_set_condition
        use filename_mod
        use path_operator_m
        double precision DTa, dt, L, U
        double precision :: direction_g(3)
        character(99) path2FlowFile
        integer i, n_unit, num_drop

        OPEN(newunit=n_unit, FILE=path%DIR//conditionFName, STATUS='OLD')
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
            read(n_unit,*) DTa
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

        Rdt = dt/DTa                       !データ読み込み時に時間軸合わせるパラメータ

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

        if(loopf - loops > 0) print*, 'Loop is from', loops, 'to', loopf
        print*, 'Delta_Time =', dt
        print*, 'Rdt', Rdt

        call set_dir_from_path(path2FlowFile, PATH_FlowDIR, FNAME_FMT)

        call check_FILE_GRID    !気流ファイルのタイプをチェック

        call set_basical_variables(dt, L, U)

        call set_gravity_acceleration(direction_g)

    end subroutine read_and_set_condition

    subroutine set_initialDroplet
        implicit none
        character(99) fname
        type(virusDroplet_t), allocatable :: droplets_ini(:)

        if(num_restart > 0) then

            print*, 'RESTRAT'

            write(fname,'("'//trim(path%backup)//'backup", i8.8, ".bu")') num_restart
            droplets = read_backup(fname)   !ここで自動割り付け

            n_start = num_restart
            n_time = n_start

        else
            droplets_ini = get_initialState_of_droplets()
            allocate(droplets(size(droplets_ini)))
            droplets(:)%virusDroplet_t = droplets_ini(:)
            n_start = 0
            n_time = n_start
            call output_droplet  !リスタートでないなら初期配置出力

        end if

        num_droplets = size(droplets)
        print*, 'num_droplets =', num_droplets
            
    end subroutine set_initialDroplet

    subroutine first_refCELLserch
        integer i, j, num_drop

        num_drop = size(droplets)

        if(unstructuredGrid) then
            do i = 1, size(leaderID) - 1
                j = leaderID(i)
                call search_refCELL(real(droplets(j)%position(:)), droplets(j)%refCELL%ID)
                droplets(leaderID(i)+1 : leaderID(i+1)-1)%refCELL%ID = droplets(j)%refCELL%ID !時間短縮を図る

                do j = leaderID(i) + 1, leaderID(i+1) - 1
                    call search_refCELL(real(droplets(j)%position(:)), droplets(j)%refCELL%ID)
                end do
            end do

        else
            do i = 1, size(leaderID) - 1
                j = leaderID(i)
                call search_refCELL_onCUBE(real(droplets(j)%position(:)), droplets(j)%refCELL)
                droplets(leaderID(i)+1 : leaderID(i+1)-1)%refCELL = droplets(j)%refCELL !時間短縮を図る

                do j = leaderID(i) + 1, leaderID(i+1) - 1
                    call search_refCELL_onCUBE(real(droplets(j)%position(:)), droplets(j)%refCELL)
                end do
            end do

        end if

    end subroutine first_refCELLserch

    subroutine survival_check
        integer vfloat, vn
        ! double precision rand
        ! double precision, save :: death_rate = 0.d0
            
        vfloat = count(droplets(:)%status == 0)
        if(vfloat == 0) then
            print*, 'No Droplet is Floating', n_time
            return  !浮遊数がゼロならリターン
        end if

        do vn = 1, num_droplets
            if ((droplets(vn)%status == 0).and.(droplets(vn)%death_param > survival_rate(n_time))) then
                droplets(vn)%status = -1
                droplets(vn)%velocity(:) = 0.0d0
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
            do vn = 1, num_droplets
                if(droplets(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call evaporation(vn)    !蒸発方程式関連の処理
                call motion_calc(vn)     !運動方程式関連の処理
            end do

        else
            do vn = 1, num_droplets
                if(droplets(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call evaporation(vn)    !蒸発方程式関連の処理
                call motion_calc_onCUBE(vn)     !運動方程式関連の処理
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

    subroutine motion_calc(vn)
        integer, intent(in) :: vn
        double precision  :: X(3), velAir(3)
        integer RefCellID
        logical stopflag, first
    
        X(:) = droplets(vn)%position(:)

        RefCellID = droplets(vn)%refCELL%ID

        stopflag = adhesion_check(vn, RefCellID)

        call area_check(X(:), stopflag)

        velAir(:) = CELLs(RefCellID)%flowVelocity(:)
    
        call calc_motionResult(droplets(vn)%virusDroplet_t, velAir, stopflag)

        call search_refCELL(real(X(:)), droplets(vn)%refCELL%ID)
        
    end subroutine motion_calc

    subroutine motion_calc_onCUBE(vn)
        use adhesion_onSTL_m
        integer, intent(in) :: vn
        double precision  :: X(3), velAir(3)
        type(reference_cell_t) :: RefC
        logical stopflag, first
    
        X(:) = droplets(vn)%position(:)

        RefC = droplets(vn)%refCELL

        stopflag = adhesion_check_onSTL(real(droplets(vn)%position(:)))

        call area_check(X(:), stopflag)

        velAir(:) = get_velocity_f(RefC%nodeID, RefC%ID)

        call calc_motionResult(droplets(vn)%virusDroplet_t, velAir, stopflag)

        call search_refCELL_onCUBE(real(X(:)), droplets(vn)%refCELL)
    
    end subroutine motion_calc_onCUBE

    subroutine calc_motionResult(droplet, Va, stopflag)
        type(virusDroplet_t), intent(inout) :: droplet
        double precision, intent(in) :: Va(3)
        logical, intent(in) :: stopflag

        if (stopflag) then
            droplet%status = 1
            droplet%velocity(:) = 0.0d0     !速度をゼロに
    
        else
            call solve_motionEquation(droplet%position(:), droplet%velocity(:), Va(:), droplet%radius)
            
        end if

    end subroutine calc_motionResult
                    
    logical function adhesion_check(vn, NCN)
        use vector_m
        integer JJ, JB
        integer, intent(in) :: vn, NCN
        double precision :: r_vector(3), inner

        adhesion_check = .false.

        do JJ = 1, size(CELLs(NCN)%boundFaceID)
            JB = CELLs(NCN)%boundFaceID(JJ)

            r_vector(:) = droplets(vn)%position(:) - BoundFACEs(JB)%center(:)

            inner = dot_product(real(r_vector(:)), BoundFACEs(JB)%normalVector(:))

            if (inner > 0.0d0) then
                adhesion_check = .true. !外向き法線ベクトルと位置ベクトルの内積が正なら付着判定
                droplets(vn)%adhesBoundID = JB     !付着した境界面番号
            end if
        end do

    end function adhesion_check

    integer function get_flowStep()
        integer Lamda, Delta

        get_flowStep = int(dble(n_time)*RDT) + OFFSET   !気流計算における経過ステップ数に相当

        Lamda = LoopF - LoopS
        
        if((Lamda > 0).and.(get_flowStep > LoopF)) then
            Delta = mod(get_flowStep - LoopS, Lamda)
            get_flowStep = LoopS + Delta
        end if

    end function

    subroutine update_flow_check
        double precision Step_air

        if(INTERVAL_FLOW <= 0) return

        Step_air = dble(n_time)*Rdt          !気流計算における経過ステップ数に相当
        if(mod(Step_air, dble(INTERVAL_FLOW)) == 0.d0) then
            call read_flow_field(first=.false.)   !流れ場の更新
        end if

    end subroutine update_flow_check

    subroutine read_flow_field(first)
        integer flowStep, FNUM
        logical, intent(in) :: first

        if(INTERVAL_FLOW <= 0) then
            call read_steadyFlowData
        else
            flowStep = get_flowStep()
            FNUM = get_FileNumber(flowStep)
            call read_unsteadyFlowData(FNUM)

        end if

        call set_MinMaxCDN

        if(first) then
            call preprocess_onFlowField         !流れ場の前処理
            if(n_start == 0) call first_refCELLserch
        end if

        if(unstructuredGrid) then
            call boundary_setting(first)
            if(.not.first) call boundary_move
        end if
            
    end subroutine read_flow_field
                      
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
                call area_check(droplets(vn)%position(:))
    
            end if
        
        end do

        ! print*, 'FIN:boundary_move'

    end subroutine boundary_move

    subroutine area_check(x, check)
        double precision, intent(inout) :: x(3)
        logical,optional,intent(inout) :: check
        integer L

        do L = 1, 3
    
            if(x(L) < MIN_CDN(L)) then
                x(L) = MIN_CDN(L)
                if(present(check)) check = .true.
            else if(x(L) > MAX_CDN(L)) then
                x(L) = MAX_CDN(L)
                if(present(check)) check = .true.
            end if

        end do

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
                drop_counter = 0

        end select

    end function drop_counter

    subroutine coalescence_check
        integer d1, d2
        integer, save :: last_coalescence = 0
        double precision :: distance, r1, r2

        if(last_coalescence == 0) last_coalescence = n_time
        !最後の合体から100ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((n_time - last_coalescence) > 100) return

        print*, 'Coalescence_check', n_time

        drop1 : do d1 = 1, num_droplets - 1
            if(droplets(d1)%status/=0) cycle drop1

            drop2 : do d2 = d1 + 1, num_droplets
                if(droplets(d2)%status/=0) cycle drop2

                distance = norm2(droplets(d2)%position(:) - droplets(d1)%position(:))
                r1 = droplets(d1)%radius
                r2 = droplets(d2)%radius

                if((r1+r2) >= distance) then
                    print*, 'Coalescence', d1, d2
                    if(r1 >= r2) then
                        call coalescence(droplets(d1)%virusDroplet_t, droplets(d2)%virusDroplet_t)
                    else
                        call coalescence(droplets(d2)%virusDroplet_t, droplets(d1)%virusDroplet_t)
                    end if
                    last_coalescence = n_time

                end if

            end do drop2

        end do drop1

    end subroutine coalescence_check

    subroutine coalescence(droplet1, droplet2)
        type(virusDroplet_t), intent(inout) :: droplet1, droplet2
        double precision :: radius_c, volume1, volume2, velocity_c(3)

        volume1 = droplet1%radius**3
        volume2 = droplet2%radius**3
        radius_c = (volume1 + volume2)**(1.d0/3.d0) !結合後の飛沫半径
        velocity_c(:) = (volume1*droplet1%velocity(:) + volume2*droplet2%velocity(:)) / (volume1 + volume2)
        
        droplet1%radius = radius_c
        droplet1%velocity(:) = velocity_c(:)
        droplet2%radius = 0.d0
        droplet2%velocity(:) = 0.d0
        droplet2%status = -2
        
    end subroutine coalescence

    subroutine output_droplet
        character(99) fname
        character(4) :: head_out = 'drop'

        write(fname,'("'//path%VTK//trim(head_out)//'",i8.8,".vtk")') n_time
        call output_droplet_VTK(fname, droplets(:)%virusDroplet_t)

        fname = path%DIR//'/particle.csv'
        call output_droplet_CSV(fname, droplets(:)%virusDroplet_t, n_time)

        call output_backup
    end subroutine

    function read_backup(fname) result(droplets_read)
        implicit none
        character(*), intent(in) :: fname
        type(Droplet_onFlow), allocatable :: droplets_read(:)
        integer i, n_unit, num_drop
    
        print*, 'READ:', fname
        open(newunit=n_unit, form='unformatted', file=fname, status='old')
            read(n_unit) num_drop

            allocate(droplets_read(num_drop))

            do i = 1, num_drop
                read(n_unit) droplets_read(i)
            end do
        close(n_unit)
    
      
    end function read_backup

    subroutine output_backup
        implicit none
        integer i, n_unit, num_drop
        character(99) fname

        num_drop = size(droplets(:))
        write(fname,'("'//path%backup//'backup", i8.8, ".bu")') n_time

        open(newunit=n_unit, form='unformatted', file=fname, status='replace')
            write(n_unit) num_drop
            do i = 1, num_drop
                write(n_unit) droplets(i)
            end do
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine output_backup

    real function environment(name)
        character(*), intent(in) :: name

        select case(name)
            case('Temperature')
                environment = T
            case('Relative Humidity')
                environment = RH
            case default
                environment = 0.0
        end select

    end function environment

    subroutine deallocation_droplet

        deallocate(droplets)
        
    end subroutine deallocation_droplet

end module drop_motion_mod