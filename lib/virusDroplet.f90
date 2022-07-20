module virusDroplet_m
    implicit none
    private

    type, public :: virusDroplet_t
        double precision :: position(3), velocity(3)=0.d0
        double precision, private :: radius, radius_min, initialRadius, deadline
        integer, private :: status=0
        integer :: coalesID=0, refCellID=0, adhesBoundID=0

        contains

        procedure :: isFloating => isDropletFloating
        procedure :: coalescenceID => dropletCoalescneceID
        procedure stop_droplet, isEvaporating, evaporation, get_radius

    end type

    type, public :: DropletGroup
        type(virusDroplet_t), allocatable :: droplet(:)

        integer, allocatable :: statusCSV(:)

        contains

        procedure output_backup
        procedure :: output_VTK => output_droplet_VTK
        procedure :: output_CSV => output_droplet_CSV

        procedure :: counter => dropletCounter
        procedure :: IDinBox => dropletIDinBox
        procedure :: inBox => dropletInBox
        procedure :: totalVolume => dropletTotalVolume
        procedure :: IDinState => dropletIDinState
        procedure :: getArea => get_dropletGroupArea

        procedure set_initialRadius, set_radiusLowerLimit
        procedure :: set_status => set_dropletGroupStatus
        procedure :: set_deadline => set_virusDeadline

        procedure survival_check
        procedure coalescence_check
        ! procedure :: calculation => Calculation_Droplets

        ! procedure :: append => append_dropletGroup

    end type

    public read_backup

    contains

    logical function isDropletFloating(self)
        class(virusDroplet_t), intent(in) :: self

        if(self%status==0) then
            isDropletFloating = .true.
        else
            isDropletFloating = .false.
        end if

    end function


    double precision function get_radius(self)
        class(virusDroplet_t), intent(in) :: self

        get_radius = self%radius

    end function

    logical function isEvaporating(self)
        class(virusDroplet_t), intent(in) :: self

        if(self%radius > self%radius_min) then
            isEvaporating = .true.
        else
            isEvaporating = .false.
        end if

    end function

    subroutine evaporation(self, dr)
        class(virusDroplet_t) self
        double precision, intent(in) :: dr

        self%radius = max(self%radius + dr, self%radius_min)

    end subroutine

    integer function dropletCoalescneceID(self)
        class(virusDroplet_t), intent(in) :: self

        dropletCoalescneceID = self%coalesID

    end function

    integer function get_statusNumber(name)
        character(*), intent(in) :: name

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

    subroutine survival_check(self, time)
        ! use terminalControler_m
        class(DropletGroup) self
        double precision, intent(in) :: time
        integer vfloat, vn
        ! double precision rand
        ! double precision, save :: death_rate = 0.d0
            
        vfloat = count(self%droplet(:)%status == 0)
        if(vfloat == 0) return  !浮遊数がゼロならリターン

        do vn = 1, size(self%droplet)
            if ((self%droplet(vn)%status == 0).and.&
                (time > self%droplet(vn)%deadline)) then

                call self%droplet(vn)%stop_droplet(status=-1)

            end if
        end do

        ! death_rate = death_rate + dble(vfloat)*(survival_rate(timeStep) - survival_rate(timeStep+1))    !このステップで死滅すべき飛沫数
        
        ! do while(death_rate >= 1.0d0)
        !     call random_number(rand)    !死滅IDは乱数で決まる
        !     vn = int(num_droplets*rand)
        !     if(vn < 1) cycle
        !     if (self%droplet(vn)%status == 0) then !浮遊粒子からのみ除去する
        !         self%droplet(vn)%status = -1
        !         ! self%droplet(vn)%position(:) = MIN_CDN(:) - 1.0d0 !計算エリア外に配置（不要かも）
        !         self%droplet(vn)%velocity(:) = 0.0d0
        !         death_rate = death_rate - 1.0d0

        !     end if
        ! end do
      
    end subroutine

    integer function dropletCounter(self, name)
        class(DropletGroup) self
        character(*), intent(in) :: name

        select case(name)
            case('total')
                dropletCounter = size(self%droplet)

            case default
                dropletCounter = count(self%droplet(:)%status == get_statusNumber(name))

        end select

    end function

    function dropletIDinBox(self, min_cdn, max_cdn, status) result(ID_array)
        class(DropletGroup) self
        double precision, intent(in) :: min_cdn(3), max_cdn(3)
        integer, intent(in), optional :: status
        double precision position(3)
        integer, allocatable :: ID_array(:)
        integer i, id_array_(size(self%droplet)), cnt

        cnt = 0
        if(present(status)) then
            do i = 1, size(self%droplet)
                if(self%droplet(i)%status /= status) cycle

                position(:) = self%droplet(i)%position(:)

                if(      ((min_cdn(1) <= position(1)) .and. (position(1) <= max_cdn(1))) &
                    .and.((min_cdn(2) <= position(2)) .and. (position(2) <= max_cdn(2))) &
                    .and.((min_cdn(3) <= position(3)) .and. (position(3) <= max_cdn(3)))    ) then

                    cnt = cnt + 1
                    id_array_(cnt) = i

                end if

            end do


        else
            do i = 1, size(self%droplet)
                position(:) = self%droplet(i)%position(:)

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

    type(DropletGroup) function dropletInBox(self, min_cdn, max_cdn)
        class(DropletGroup) self
        double precision, intent(in) :: min_cdn(3), max_cdn(3)
        integer, allocatable :: IDinBox(:)

        IDinBox = self%IDinBox(min_cdn, max_cdn)
        dropletInBox%droplet = self%droplet(IDinBox)
        
    end function

    double precision function dropletTotalVolume(self)
        class(DropletGroup) self
        integer i
        double precision, parameter :: PI = acos(-1.d0) 

        dropletTotalVolume = 0.d0

        do i = 1, size(self%droplet)
            dropletTotalVolume = dropletTotalVolume + self%droplet(i)%initialRadius**3
        end do

        dropletTotalVolume = dropletTotalVolume * 4.d0/3.d0*PI

    end function

    function dropletIDinState(self, status) result(ID_array)
        class(DropletGroup) self
        character(*), intent(in) :: status
        integer, allocatable :: ID_array(:)
        integer i, cnt, statusNumber

        cnt = 0
        statusNumber = get_statusNumber(status)
        allocate(ID_array(count(self%droplet(:)%status==statusNumber)))
        
        do i = 1, size(self%droplet)
            if(self%droplet(i)%status == statusNumber) then

                cnt = cnt + 1
                ID_array(cnt) = i

            end if

        end do
        
    end function

    subroutine get_dropletGroupArea(self, AreaMin, AreaMax)
        class(DropletGroup) self
        double precision, intent(out) :: AreaMin(3), AreaMax(3)
        integer i, m

        AreaMin(:) = 1.d9
        AreaMax(:) = -1.d9
        do m = 1, size(self%droplet)
            if(self%droplet(m)%status==0) then
                do i = 1, 3
                    AreaMin(i) = min(self%droplet(m)%position(i), AreaMin(i))
                    AreaMax(i) = max(self%droplet(m)%position(i), AreaMax(i))
                end do
            end if
        end do

    end subroutine

    subroutine coalescence_check(self, stat)
        class(DropletGroup) self
        integer, intent(out), optional :: stat
        integer d1, d2, num_droplets, num_coales
        double precision :: distance, r1, r2

        num_coales = 0

        num_droplets = size(self%droplet)

        !$OMP parallel do private(distance, r1, r2)
        drop1 : do d1 = 1, num_droplets - 1
            if(self%droplet(d1)%status/=0) cycle drop1

            r1 = self%droplet(d1)%radius

            drop2 : do d2 = d1 + 1, num_droplets
                if(self%droplet(d2)%status/=0) cycle drop2

                r2 = self%droplet(d2)%radius

                distance = norm2(self%droplet(d2)%position(:) - self%droplet(d1)%position(:))

                if((r1+r2) >= distance) then
                    ! print*, d1, 'and', d2, 'coalesce!'
                    if(r1 >= r2) then
                        call coalescence(self%droplet(d1), self%droplet(d2), d1)
                    else
                        call coalescence(self%droplet(d2), self%droplet(d1), d2)
                    end if
                    num_coales = num_coales + 1

                end if

            end do drop2

        end do drop1
        !$OMP end parallel do

        if(present(stat)) stat = num_coales

    end subroutine coalescence_check

    subroutine coalescence(droplet1, droplet2, baseID)
        type(virusDroplet_t), intent(inout) :: droplet1, droplet2
        integer, intent(in) :: baseID
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
        droplet2%status = -2
        droplet2%coalesID = baseID
        
    end subroutine coalescence

    subroutine output_backup(self, fname)
        class(DropletGroup) self
        character(*), intent(in) :: fname
        integer i, n_unit, num_drop

        num_drop = size(self%droplet(:))

        open(newunit=n_unit, form='unformatted', file=fname, status='replace')
            write(n_unit) num_drop
            do i = 1, num_drop
                write(n_unit) self%droplet(i)
            end do
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine

    subroutine output_droplet_VTK(self, fname, deadline)
        class(DropletGroup) self
        character(*), intent(in) :: fname
        logical, optional :: deadline
        integer vn, n_unit, num_drop

        num_drop = size(self%droplet)
        
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I0,A)') 'POINTS ', num_drop,' float'                              !節点の数
            DO vn = 1, num_drop                                                            !節点の数だけループ
                write(n_unit,'(3(f10.5,2X))') self%droplet(vn)%position(:)   !節点の座標（左から順にx,y,z）
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
                write(n_unit,'(I0)') self%droplet(vn)%status                                            !飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO

            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(e10.3)') self%droplet(vn)%radius*2.0d0                            !飛沫の直径
            END DO

            if(present(deadline)) then
                if(deadline) then
                    write(n_unit,'(A)') 'SCALARS Deadline float'
                    write(n_unit,'(A)') 'LOOKUP_TABLE default'
                    DO vn = 1, num_drop                                                            !セルの数だけループ
                        write(n_unit,'(e10.3)') self%droplet(vn)%deadline
                    END DO
                end if
            end if

            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1, num_drop                                                             !セルの数だけループ
                write(n_unit,'(3(f10.5,2X))') self%droplet(vn)%velocity(:)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine output_droplet_VTK

    subroutine output_droplet_CSV(self, fname, time, initial)
        class(DropletGroup) self
        character(*), intent(in) :: fname
        double precision, intent(in) :: time
        logical, intent(in) :: initial
        integer n_unit, J

        if(initial) then !初回ならファイル新規作成
            open(newunit=n_unit, file=fname, status='replace')

        else
            open(newunit=n_unit, file=fname, action='write', status='old', position='append')

        end if

        if(.not.allocated(self%statusCSV)) self%statusCSV = [0, 1, -1,-2]
        
        write(n_unit,'(*(g0:,","))') real(time), (count(self%droplet(:)%status==self%statusCSV(J)), J = 1, size(self%statusCSV))

        close(n_unit)

    end subroutine

    ! subroutine append_dropletGroup(self, dGroup)
    !     class(DropletGroup) self
    !     type(DropletGroup), intent(in) :: dGroup

    !     self%droplet = [self%droplet, dGroup%droplet]

    ! end subroutine

    function read_backup(fname) result(dGroup_read)
        character(*), intent(in) :: fname
        type(DropletGroup) dGroup_read
        integer i, n_unit, num_drop
    
        print*, 'READ : ', trim(fname)
        open(newunit=n_unit, form='unformatted', file=fname, status='old', action='read')
            read(n_unit) num_drop

            allocate(dGroup_read%droplet(num_drop))

            do i = 1, num_drop
                read(n_unit) dGroup_read%droplet(i)
            end do
        close(n_unit)
    
    end function
    
    function read_droplet_VTK(fname) result(dGroup_read)
        implicit none
        character(*), intent(in) :: fname
        type(DropletGroup) dGroup_read
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

            allocate(dGroup_read%droplet(num_drop), diameter(num_drop))

            DO vn = 1, num_drop
                read(n_unit, *) dGroup_read%droplet(vn)%position(:)
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
                read(n_unit, *) dGroup_read%droplet(vn)%status
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) diameter(vn)
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) dGroup_read%droplet(vn)%deadline
            END DO

            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) dGroup_read%droplet(vn)%velocity(:)
            END DO
        close(n_unit)
    
        dGroup_read%droplet(:)%radius = diameter(:) * 0.5d0
      
    end function

    subroutine stop_droplet(self, status)
        class(virusDroplet_t) self
        integer, optional :: status

        self%velocity(:) = 0.0d0
        if(present(status)) then
            self%status = status
        else
            self%status = 1
        end if

    end subroutine

    subroutine set_initialRadius(self, radius)
        class(DropletGroup) self
        double precision, intent(in) :: radius(:)

        self%droplet(:)%initialRadius = radius(:)
        self%droplet(:)%radius = self%droplet(:)%initialRadius

    end subroutine

    subroutine set_radiusLowerLimit(self, lowerLimitRatio)
        class(DropletGroup) self
        double precision, intent(in) :: lowerLimitRatio

        self%droplet(:)%radius_min = self%droplet(:)%initialRadius*lowerLimitRatio

    end subroutine

    subroutine set_virusDeadline(self, deadline)
        class(DropletGroup) self
        double precision, intent(in) :: deadline(:)

        self%droplet(:)%deadline = deadline(:)

    end subroutine

    subroutine set_dropletGroupStatus(self, status, ID)
        class(DropletGroup) self
        character(*), intent(in) :: status
        integer, intent(in), optional :: ID(:)

        if(present(ID)) then
            self%droplet(ID)%status = get_statusNumber(status)
        else
            self%droplet(:)%status = get_statusNumber(status)
        end if

    end subroutine

end module virusDroplet_m