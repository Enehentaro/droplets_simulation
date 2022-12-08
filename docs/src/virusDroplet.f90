module virusDroplet_m
    implicit none
    private

    integer, parameter :: statusCSV(4) = [0, 1, -1,-2]

    type, public :: virusDroplet_t
        !! ウイルス飛沫構造体

        double precision :: position(3), velocity(3)=0.d0
        double precision, private :: radius, radius_min, initialRadius, deadline
        integer, private :: status=0
        integer :: coalesID=0, refCellID=0, adhesBoundID=0

        contains

        procedure :: isFloating => isDropletFloating
        procedure :: coalescenceID => dropletCoalescneceID
        procedure stop_droplet, isEvaporating, evaporation, get_radius

    end type

    public read_backup
    public dropletCounter, dropletIDinState, dropletIDinBox
    public set_dropletStatus, set_initialRadius, set_radiusLowerLimit, set_virusDeadline
    public survival_check
    public output_backup, output_droplet_CSV, output_droplet_VTK
    public get_dropletsArea, dropletTotalVolume
    public coalescence_check

    contains

    logical function isDropletFloating(self)
        !! 飛沫が浮遊しているか否かを返す
        class(virusDroplet_t), intent(in) :: self

        if(self%status==0) then
            isDropletFloating = .true.
        else
            isDropletFloating = .false.
        end if

    end function


    double precision function get_radius(self)
        !! 現時刻の飛沫半径を返す
        class(virusDroplet_t), intent(in) :: self

        get_radius = self%radius

    end function

    logical function isEvaporating(self)
        !! 飛沫が限界まで蒸発したか否かを返す
        class(virusDroplet_t), intent(in) :: self

        if(self%radius > self%radius_min) then
            isEvaporating = .true.
        else
            isEvaporating = .false.
        end if

    end function

    subroutine evaporation(self, dr)
        !! 蒸発の処理を行う
        class(virusDroplet_t) self
        double precision, intent(in) :: dr
            !!半径変化量（蒸発ならマイナスの値を指定）

        self%radius = max(self%radius + dr, self%radius_min)

    end subroutine

    integer function dropletCoalescneceID(self)
        !! 合体飛沫の合体先のIDを返す
        class(virusDroplet_t), intent(in) :: self

        dropletCoalescneceID = self%coalesID

    end function

    integer function get_statusNumber(name)
        !! 飛沫の状態を対応する整数値で返す
        character(*), intent(in) :: name
            !! 状態

        select case(name)
            case('floating')
                get_statusNumber = 0

            case('adhesion')
                get_statusNumber = 1

            case('death')
                get_statusNumber = -1

            case('coalescence')
                get_statusNumber = -2

            case('nonActive')
                get_statusNumber = -99

            case default
                print*, '**ERROR [statusNumber] : ', name, ' is not found.**'
                error stop

        end select

    end function

    subroutine survival_check(droplets, time)
        !! 飛沫の不活性化判定
        ! use terminalControler_m
        type(virusDroplet_t) droplets(:)
        double precision, intent(in) :: time
            !! 現在時刻（無次元）
        integer vn
        ! double precision rand
        ! double precision, save :: death_rate = 0.d0

        do vn = 1, size(droplets)
            !浮遊中であり、時刻が寿命以上であれば不活性化
            if ((droplets(vn)%isFloating()) .and. (time > droplets(vn)%deadline)) then

                call droplets(vn)%stop_droplet(status='death')

            end if
        end do

        ! death_rate = death_rate + dble(vfloat)*(survival_rate(timeStep) - survival_rate(timeStep+1))    !このステップで死滅すべき飛沫数
        
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
      
    end subroutine

    integer function dropletCounter(droplets, name)
        !! 指定された状態にある飛沫数をカウント
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: name
            !! キーワード（total, floating, etc.）

        select case(name)
            case('total')
                dropletCounter = size(droplets)

            case default
                dropletCounter = count(droplets(:)%status == get_statusNumber(name))

        end select

    end function

    function dropletIDinBox(droplets, min_cdn, max_cdn, status) result(ID_array)
        !! 直方体領域内の飛沫のID配列を返す
        type(virusDroplet_t) droplets(:)

        double precision, intent(in) :: min_cdn(3)
            !! 直方体領域の最小座標
        double precision, intent(in) :: max_cdn(3)
            !! 直方体領域の最大座標

        integer, intent(in), optional :: status
            !! 状態

        double precision position(3)
        integer, allocatable :: ID_array(:)
        integer i, id_array_(size(droplets)), cnt

        cnt = 0
        if(present(status)) then
            do i = 1, size(droplets)
                if(droplets(i)%status /= status) cycle

                position(:) = droplets(i)%position(:)

                if(      ((min_cdn(1) <= position(1)) .and. (position(1) <= max_cdn(1))) &
                    .and.((min_cdn(2) <= position(2)) .and. (position(2) <= max_cdn(2))) &
                    .and.((min_cdn(3) <= position(3)) .and. (position(3) <= max_cdn(3)))    ) then

                    cnt = cnt + 1
                    id_array_(cnt) = i

                end if

            end do


        else
            do i = 1, size(droplets)
                position(:) = droplets(i)%position(:)

                if(      ((min_cdn(1) <= position(1)) .and. (position(1) <= max_cdn(1))) &
                    .and.((min_cdn(2) <= position(2)) .and. (position(2) <= max_cdn(2))) &
                    .and.((min_cdn(3) <= position(3)) .and. (position(3) <= max_cdn(3)))    ) then

                    cnt = cnt + 1
                    id_array_(cnt) = i

                end if

            end do

        end if

        ID_array = id_array_(:cnt)
        
    end function

    double precision function dropletTotalVolume(droplets)
        !! 配列内の全飛沫の総体積（無次元）を計算

        type(virusDroplet_t) droplets(:)
        integer i
        double precision, parameter :: PI = acos(-1.d0) 

        dropletTotalVolume = 0.d0

        do i = 1, size(droplets)
            dropletTotalVolume = dropletTotalVolume + droplets(i)%initialRadius**3
        end do

        dropletTotalVolume = dropletTotalVolume * 4.d0/3.d0*PI

    end function

    function dropletIDinState(droplets, status) result(ID_array)
        !! 任意の状態の飛沫のID配列を返す
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: status
        integer, allocatable :: ID_array(:)
        integer i, cnt, statusNumber

        cnt = 0
        statusNumber = get_statusNumber(status)
        allocate(ID_array(count(droplets(:)%status==statusNumber)))
        
        do i = 1, size(droplets)
            if(droplets(i)%status == statusNumber) then

                cnt = cnt + 1
                ID_array(cnt) = i

            end if

        end do
        
    end function

    subroutine get_dropletsArea(droplets, AreaMin, AreaMax)
        !! 飛沫配列内の飛沫の座標の最大最小を返す
        type(virusDroplet_t) droplets(:)
        double precision, intent(out) :: AreaMin(3), AreaMax(3)
        integer i, m

        AreaMin(:) = 1.d9
        AreaMax(:) = -1.d9
        do m = 1, size(droplets)
            if(droplets(m)%status==0) then
                do i = 1, 3
                    AreaMin(i) = min(droplets(m)%position(i), AreaMin(i))
                    AreaMax(i) = max(droplets(m)%position(i), AreaMax(i))
                end do
            end if
        end do

    end subroutine

    subroutine coalescence_check(droplets, stat)
        !! 合体判定
        type(virusDroplet_t) droplets(:)
        integer, intent(out), optional :: stat
        integer d1, d2, num_droplets, num_coales
        double precision :: distance, r1, r2

        num_coales = 0

        num_droplets = size(droplets)

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
                        call coalescence(droplets(d1), droplets(d2), d1)
                    else
                        call coalescence(droplets(d2), droplets(d1), d2)
                    end if
                    num_coales = num_coales + 1

                end if

            end do drop2

        end do drop1
        !$OMP end parallel do

        if(present(stat)) stat = num_coales

    end subroutine

    subroutine coalescence(droplet1, droplet2, baseID)
        !! 合体処理
        type(virusDroplet_t), intent(inout) :: droplet1, droplet2
        integer, intent(in) :: baseID
            !! 合体先（ベースとなる飛沫）のID
        double precision r3_1, r3_2, position_c(3), velocity_c(3)

        r3_1 = droplet1%radius**3
        r3_2 = droplet2%radius**3
        position_c(:) = (r3_1*droplet1%position(:) + r3_2*droplet2%position(:)) / (r3_1 + r3_2)
        velocity_c(:) = (r3_1*droplet1%velocity(:) + r3_2*droplet2%velocity(:)) / (r3_1 + r3_2)
        
        droplet1%radius = (r3_1 + r3_2)**(1.d0/3.d0)
        droplet1%position(:) = position_c(:)
        droplet1%velocity(:) = velocity_c(:)
        droplet1%radius_min = (droplet1%radius_min**3 + droplet2%radius_min**3)**(1.d0/3.d0)
        ! droplet1%initialRadius = radius_afterCoalescence(droplet1%initialRadius, droplet2%initialRadius)

        droplet2%radius = 0.d0
        droplet2%position(:) = position_c(:)
        droplet2%velocity(:) = velocity_c(:)
        droplet2%status = get_statusNumber('coalescence')
        droplet2%coalesID = baseID
        
    end subroutine

    subroutine output_backup(droplets, fname)
        !! backupファイルの出力。
        !! 配列をループ使わずそのまま書き出すほうが多分ファイルサイズ効率が良いので、いつか修正したい。
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: fname
        integer i, n_unit, num_drop

        num_drop = size(droplets(:))

        open(newunit=n_unit, form='unformatted', file=fname, status='replace')
            write(n_unit) num_drop
            do i = 1, num_drop
                write(n_unit) droplets(i)
            end do
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine

    subroutine output_droplet_VTK(droplets, fname, deadline)
        !! VTK形式でファイル出力
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: fname
        logical, optional :: deadline
        integer vn, n_unit, num_drop

        num_drop = size(droplets)
        
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I0,A)') 'POINTS ', num_drop,' float'                              !節点の数
            DO vn = 1, num_drop                                                            !節点の数だけループ
                write(n_unit,'(3(f10.5,2X))') droplets(vn)%position(:)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I0,2X,I0)') 'CELLS ', num_drop, num_drop*2                          !セルの数、セルの数×2
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I0,2X,I0)') 1, vn-1                                          !セルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I0)') 'CELL_TYPES ', num_drop
            DO vn = 1, num_drop                                                           !セルの数だけループ
                write(n_unit,'(I0)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I0)') 'CELL_DATA ', num_drop                                      !ここからセルのデータという合図、セルの数
            
            write(n_unit,'(A)') 'SCALARS Status int'                                   !飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I0)') droplets(vn)%status                                            !飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO

            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(e10.3)') droplets(vn)%radius*2.0d0                            !飛沫の直径
            END DO

            if(present(deadline)) then
                if(deadline) then
                    write(n_unit,'(A)') 'SCALARS Deadline float'
                    write(n_unit,'(A)') 'LOOKUP_TABLE default'
                    DO vn = 1, num_drop                                                            !セルの数だけループ
                        write(n_unit,'(e10.3)') droplets(vn)%deadline
                    END DO
                end if
            end if

            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1, num_drop                                                             !セルの数だけループ
                write(n_unit,'(3(f10.5,2X))') droplets(vn)%velocity(:)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine

    subroutine output_droplet_CSV(droplets, fname, time, initial)
        !! CSV形式で時系列データ出力
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: fname
        double precision, intent(in) :: time
        logical, intent(in) :: initial
        integer n_unit, J

        if(initial) then !初回ならファイル新規作成
            open(newunit=n_unit, file=fname, status='replace')

        else
            open(newunit=n_unit, file=fname, action='write', status='old', position='append')

        end if
        
        write(n_unit,'(*(g0:,","))') real(time), (count(droplets(:)%status==statusCSV(J)), J = 1, size(statusCSV))

        close(n_unit)

    end subroutine

    ! subroutine append_dropletGroup(droplets, dGroup)
    !     type(virusDroplet_t) droplets(:)
    !     type(DropletGroup), intent(in) :: dGroup

    !     droplets = [droplets, dGroup%droplet]

    ! end subroutine

    function read_backup(fname) result(droplets)
        !! backupファイルを読み込み、飛沫配列を返す
        character(*), intent(in) :: fname
        type(virusDroplet_t), allocatable :: droplets(:)
        integer i, n_unit, num_drop
    
        print*, 'READ : ', trim(fname)
        open(newunit=n_unit, form='unformatted', file=fname, status='old', action='read')
            read(n_unit) num_drop

            allocate(droplets(num_drop))

            do i = 1, num_drop
                read(n_unit) droplets(i)
            end do
        close(n_unit)
    
    end function
    
    function read_droplet_VTK(fname) result(droplets)
        !! VTKファイルを読み込み、飛沫配列を返す
        implicit none
        character(*), intent(in) :: fname
        type(virusDroplet_t), allocatable :: droplets(:)
        double precision, allocatable :: diameter(:)
        integer vn, n_unit, num_drop
        character(10) str
    
        print*, 'READ : ', fname
        open(newunit=n_unit, file=fname, status='old', action='read')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit, *) str, num_drop

            allocate(droplets(num_drop), diameter(num_drop))

            DO vn = 1, num_drop
                read(n_unit, *) droplets(vn)%position(:)
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit,'()')
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit,'()')
            END DO
            read(n_unit,'()')

            read(n_unit,'()')   !CELL_DATA

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets(vn)%status
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) diameter(vn)
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets(vn)%deadline
            END DO

            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets(vn)%velocity(:)
            END DO
        close(n_unit)
    
        droplets(:)%radius = diameter(:) * 0.5d0
      
    end function

    subroutine stop_droplet(self, status)
        !! 飛沫を静止させる
        class(virusDroplet_t) self
        character(*), intent(in), optional :: status

        self%velocity(:) = 0.0d0
        if(present(status)) then
            self%status = get_statusNumber(status)
        else
            self%status = 1
        end if

    end subroutine

    subroutine set_initialRadius(droplets, radius)
        !! 飛沫の初期半径を配列からセットする
        type(virusDroplet_t) droplets(:)
        double precision, intent(in) :: radius(:)

        droplets(:)%initialRadius = radius(:)
        droplets(:)%radius = droplets(:)%initialRadius

    end subroutine

    subroutine set_radiusLowerLimit(droplets, lowerLimitRatio)
        !! 飛沫の蒸発限界半径をセットする
        type(virusDroplet_t) droplets(:)
        double precision, intent(in) :: lowerLimitRatio
            !! 限界半径と初期半径の比

        droplets(:)%radius_min = droplets(:)%initialRadius*lowerLimitRatio

    end subroutine

    subroutine set_virusDeadline(droplets, deadline)
        !! 飛沫の寿命を配列からセットする
        type(virusDroplet_t) droplets(:)
        double precision, intent(in) :: deadline(:)

        droplets(:)%deadline = deadline(:)

    end subroutine

    subroutine set_dropletStatus(droplets, status, ID)
        !! 飛沫の状態を一気にセットする
        !! 特定のIDの飛沫だけセットしたい場合は、ID配列を引数に渡す
        type(virusDroplet_t) droplets(:)
        character(*), intent(in) :: status
        integer, intent(in), optional :: ID(:)

        if(present(ID)) then
            droplets(ID)%status = get_statusNumber(status)
        else
            droplets(:)%status = get_statusNumber(status)
        end if

    end subroutine

end module virusDroplet_m