module dropletGroup_m
    use virusDroplet_m
    implicit none

    real, private :: T, RH

    integer, target :: timeStep = 0  !時間ステップ

    integer, allocatable :: statusCSV(:)

    type dropletGroup
        type(virusDroplet_t), allocatable :: droplet(:)

        contains

        procedure calc_initialPosition
        procedure calc_initialRadius
        procedure set_virusDeadline
        procedure calc_minimumRadius
        procedure first_refCellSearch

        procedure output_backup
        procedure :: output_VTK => output_droplet_VTK
        procedure :: output_CSV => output_droplet_CSV

        procedure :: counter => dropletCounter

        procedure adhesion_check
        procedure survival_check
        procedure coalescence_check
        procedure :: calculation => Calculation_Droplets
        procedure boundary_move

        procedure :: append => append_dropletGroup

    end type

    contains

    type(dropletGroup) function generate_dropletGroup(num_droplet)   !コンストラクタ
        integer, intent(in) :: num_droplet

        allocate(generate_dropletGroup%droplet(num_droplet))

        call generate_dropletGroup%calc_initialPosition()
        call generate_dropletGroup%calc_initialRadius()
        call generate_dropletGroup%set_virusDeadline()
        generate_dropletGroup%droplet(:)%radius = generate_dropletGroup%droplet(:)%initialRadius
        call generate_dropletGroup%calc_minimumRadius(RH) !最小半径の計算

        call generate_dropletGroup%first_refCellSearch()

    end function

    type(dropletGroup) function read_InitialDistribution(dir)
        character(*), intent(in) :: dir

        read_InitialDistribution = read_droplet_VTK(dir//'/InitialDistribution.vtk')
        read_InitialDistribution%droplet(:)%initialRadius = read_InitialDistribution%droplet(:)%radius
        call read_InitialDistribution%calc_minimumRadius(RH)

        call read_InitialDistribution%first_refCellSearch()

    end function

    !====================ここからメソッド====================

    subroutine calc_initialRadius(self)
        use csv_reader
        class(dropletGroup) self
        integer vn, i, num_drop
        double precision, allocatable :: threshold(:,:)
        integer, allocatable :: rad_cnt(:)
        double precision :: radius_dim(size(self%droplet))
        real(8) random_rad

        num_drop = size(self%droplet)
        call read_CSV('data/radius_distribution.csv', threshold)
        ! threshold = read_csv_dble('radius_distribution.csv')

        allocate(rad_cnt(size(threshold, dim=2)), source=0)

        do vn = 1, num_drop                       !飛沫半径の分布を乱数によって与える

            call random_number(random_rad)

            do i = 1, size(threshold, dim=2)
                if(random_rad < threshold(2, i)) then
                    radius_dim(vn) = threshold(1, i) * 1.0d-6
                    rad_cnt(i) = rad_cnt(i) + 1
                    exit
                end if
            end do

        end do

        self%droplet%initialRadius = radius_dim(:) / representative_value('Length') !初期飛沫半径のセットおよび無次元化

        if (sum(rad_cnt) /= num_drop) then
            print*, 'random_rad_ERROR', sum(rad_cnt), num_drop
            stop
        end if

        do i = 1, size(rad_cnt)
            print*,'rad_cnt(', threshold(1, i), ') =', rad_cnt(i)
        end do

    end subroutine

    subroutine calc_initialPosition(self)
        use csv_reader
        use filename_mod
        use caseNameList_m
        class(dropletGroup) self
        integer kx,ky,kz, num_perEdge, num_per_box, k, k_end, cnt
        integer i_box, num_box, num_drop
        double precision :: standard(3), delta(3), width(3), randble(3)
        double precision, allocatable :: position_mat(:,:)
        
        num_drop = size(self%droplet)

        call read_CSV(get_caseName()//'/'//IniPositionFName, position_mat)

        num_box = size(position_mat, dim=2)

        num_per_box = num_drop / num_box
        
        print*, 'calc_initialPosition'

        num_perEdge = 1
        do while(num_box*((num_perEdge+1)**3) < num_drop)
            num_perEdge = num_perEdge + 1    !配置帯一辺当たりの飛沫数
        end do

        k = 1
        cnt = 1
        do i_box = 1, num_box
  
            width(:) = position_mat(4:6, i_box)
            standard(:) = position_mat(1:3, i_box) - 0.5d0*width(:)
            delta(:) = width(:) / dble(num_perEdge - 1)

            do kx = 1, num_perEdge

                do ky = 1, num_perEdge

                    do kz = 1, num_perEdge

                        self%droplet(k)%position(1) = standard(1) + delta(1)*dble(kx - 1)
                        self%droplet(k)%position(2) = standard(2) + delta(2)*dble(ky - 1)
                        self%droplet(k)%position(3) = standard(3) + delta(3)*dble(kz - 1)
                        k = k + 1
                        
                    end do

                end do

            end do

            k_end = i_box * num_per_box
            if(i_box == num_box) k_end = num_drop

            do while(k <= k_end)
                call random_number(randble(:))
                self%droplet(k)%position(:) = standard(:) + width(:)*randble(:)
                k = k + 1
            end do

            print*, 'BOX', i_box, 'has', k - cnt, 'droplets.'

            cnt = k

        end do

    end subroutine

    subroutine set_virusDeadline(self)
        class(dropletGroup) self
        double precision randble(size(self%droplet))

        call random_number(randble)

        self%droplet(:)%deadline = virusDeadline(randble(:)) + Time_onSimulation(timeStep)

    end subroutine set_virusDeadline

    subroutine calc_minimumRadius(self, RelativeHumidity)
        class(dropletGroup) self
        real, intent(in) :: RelativeHumidity

        self%droplet(:)%radius_min = get_minimumRadius(self%droplet(:)%initialRadius, RelativeHumidity) !最小半径の計算

    end subroutine

    subroutine first_refCellSearch(self)
        class(dropletGroup) self
        integer j, num_drop
        logical success

        print*, 'first_refCellSearch occured!'

        num_drop = size(self%droplet)

        if(unstructuredGrid) then
            j = 1
            self%droplet(j)%refCELL%ID = nearest_cell(real(self%droplet(j)%position(:)))

            self%droplet(j+1:)%refCELL%ID = self%droplet(j)%refCELL%ID !時間短縮を図る

            do j = 2, num_drop
                call search_refCELL(real(self%droplet(j)%position(:)), self%droplet(j)%refCELL%ID, stat=success)
                if(.not.success) self%droplet(j+1:)%refCELL%ID = self%droplet(j)%refCELL%ID
            end do

        else
            j = 1
            self%droplet(j)%refCELL%ID = get_cube_contains(real(self%droplet(j)%position(:)))    
            self%droplet(j)%refCELL%nodeID(:) = nearest_node(real(self%droplet(j)%position(:)), self%droplet(j)%refCELL%ID)

            self%droplet(j+1:)%refCELL = self%droplet(j)%refCELL !時間短縮を図る

            do j = 2, num_drop
                call search_refCELL_onCUBE(real(self%droplet(j)%position(:)), self%droplet(j)%refCELL)
            end do

        end if

    end subroutine first_refCellSearch

    subroutine adhesion_check(self)
        use adhesion_onSTL_m
        class(dropletGroup) self
        integer i
        
        if(unstructuredGrid) then
            do i = 1, size(self%droplet)
                if(self%droplet(i)%status==0) then
                    call self%droplet(i)%adhesion_onBound()
                    call self%droplet(i)%area_check()
                end if
            end do
        else
            do i = 1, size(self%droplet)
                if(self%droplet(i)%status==0) then
                    if(adhesion_onSTL(real(self%droplet(i)%position(:)))) call stop_droplet(self%droplet(i))
                    call self%droplet(i)%area_check()
                end if
            end do
        end if

    end subroutine adhesion_check

    subroutine survival_check(self)
        use terminalControler_m
        class(dropletGroup) self
        integer vfloat, vn
        ! double precision rand
        ! double precision, save :: death_rate = 0.d0
            
        vfloat = count(self%droplet(:)%status == 0)
        if(vfloat == 0) return  !浮遊数がゼロならリターン

        do vn = 1, size(self%droplet)
            if ((self%droplet(vn)%status == 0).and.&
                (Time_onSimulation(timeStep) > self%droplet(vn)%deadline)) then

                call stop_droplet(self%droplet(vn), status=-1)

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
      
    end subroutine survival_check

    subroutine Calculation_Droplets(self)
        class(dropletGroup) self
        integer vn

        if(unstructuredGrid) then
            !$omp parallel do
            do vn = 1, size(self%droplet)
                if(self%droplet(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call self%droplet(vn)%evaporation()    !蒸発方程式関連の処理
                call self%droplet(vn)%motionCalculation()     !運動方程式関連の処理
            end do
            !$omp end parallel do

        else
            do vn = 1, size(self%droplet)
                if(self%droplet(vn)%status /= 0) cycle !浮遊状態でないなら無視
                call self%droplet(vn)%evaporation()    !蒸発方程式関連の処理
                call self%droplet(vn)%motionCalculation_onCUBE()     !運動方程式関連の処理
            end do

        end if

    end subroutine Calculation_Droplets
                      
    subroutine boundary_move(self) !境界面の移動に合わせて付着飛沫も移動
        class(dropletGroup) self
        integer vn, JB

        ! print*, 'CALL:boundary_move'

        do vn = 1, size(self%droplet)
        
            if (self%droplet(vn)%status <= 0) cycle !付着していないならスルー
    
            JB = self%droplet(vn)%adhesBoundID
            if (JB > 0) then
                self%droplet(vn)%position(:) &
                    = self%droplet(vn)%position(:) + BoundFACEs(JB)%moveVector(:) !面重心の移動量と同じだけ移動
            else
                call area_check(self%droplet(vn))
    
            end if
        
        end do

        ! print*, 'FIN:boundary_move'

    end subroutine boundary_move

    integer function dropletCounter(self, name)
        class(dropletGroup) self
        character(*), intent(in) :: name

        select case(name)
            case('total')
                dropletCounter = size(self%droplet)

            case('adhesion')
                dropletCounter = count(self%droplet(:)%status > 0)

            case('floating')
                dropletCounter = count(self%droplet(:)%status == 0)

            case('death')
                dropletCounter = count(self%droplet(:)%status == -1)

            case('coalescence')
                dropletCounter = count(self%droplet(:)%status == -2)

            case default
                print*, '**ERROR [dropletCounter] : ', name, ' is not found.**'
                stop

        end select

    end function

    subroutine coalescence_check(self, stat)
        class(dropletGroup) self
        integer, intent(out) :: stat
        integer d1, d2, num_droplets
        double precision :: distance, r1, r2

        stat = 0

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
                        call coalescence(self%droplet(d1), self%droplet(d2))
                    else
                        call coalescence(self%droplet(d2), self%droplet(d1))
                    end if
                    stat = stat + 1

                end if

            end do drop2

        end do drop1
        !$OMP end parallel do

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
        call droplet2%stop_droplet(status=-2)
        
    end subroutine coalescence

    subroutine output_backup(self, dir)
        class(dropletGroup) self
        character(*), intent(in) :: dir
        integer i, n_unit, num_drop
        character(99) fname

        num_drop = size(self%droplet(:))
        write(fname,'("'//dir//'/backup_", i0, ".bu")') timeStep

        open(newunit=n_unit, form='unformatted', file=fname, status='replace')
            write(n_unit) num_drop
            do i = 1, num_drop
                write(n_unit) self%droplet(i)
            end do
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine output_backup

    subroutine output_droplet_VTK(self, fname, initial)
        class(dropletGroup) self
        character(*), intent(in) :: fname
        logical, optional :: initial
        integer vn, n_unit, num_drop

        num_drop = size(self%droplet)
        
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ', num_drop,' float'                              !節点の数
            DO vn = 1, num_drop                                                            !節点の数だけループ
                write(n_unit,'(3(f20.15,2X))') self%droplet(vn)%position(:)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12,2X,I12)') 'CELLS ', num_drop, num_drop*2                          !セルの数、セルの数×2
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I4,2X,I12)') 1, vn-1                                          !セルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12)') 'CELL_TYPES', num_drop
            DO vn = 1, num_drop                                                           !セルの数だけループ
                write(n_unit,'(I4)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12)') 'CELL_DATA ', num_drop                                      !ここからセルのデータという合図、セルの数
            
            write(n_unit,'(A)') 'SCALARS Status int'                                   !飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I4)') self%droplet(vn)%status                                            !飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO

            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(f20.15)') self%droplet(vn)%radius*2.0d0                            !飛沫の直径
            END DO

            if(present(initial)) then
                if(initial) then
                    write(n_unit,'(A)') 'SCALARS Deadline float'
                    write(n_unit,'(A)') 'LOOKUP_TABLE default'
                    DO vn = 1, num_drop                                                            !セルの数だけループ
                        write(n_unit,'(e20.12)') self%droplet(vn)%deadline
                    END DO
                end if
            end if

            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1, num_drop                                                             !セルの数だけループ
                write(n_unit,'(3(f20.15,2X))') self%droplet(vn)%velocity(:)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

    end subroutine output_droplet_VTK

    subroutine output_droplet_CSV(self, fname, time, initial)
        class(dropletGroup) self
        character(*), intent(in) :: fname
        double precision, intent(in) :: time
        logical, intent(in) :: initial
        integer n_unit, L

        if(initial) then !初回ならファイル新規作成
            open(newunit=n_unit, file=fname, status='replace')

        else
            open(newunit=n_unit, file=fname, action='write', status='old', position='append')

        end if

        if(.not.allocated(statusCSV)) statusCSV = [0, 1, -1,-2]
        
        write(n_unit,'(*(g0:,","))') real(time), (count(self%droplet(:)%status==statusCSV(L)), L = 1, size(statusCSV))

        close(n_unit)

    end subroutine

    subroutine append_dropletGroup(self, dGroup)
        class(dropletGroup) self
        type(dropletGroup), intent(in) :: dGroup

        self%droplet = [self%droplet, dGroup%droplet]

    end subroutine

    !====================メソッドここまで====================

    subroutine set_environment(Temperature, RelativeHumidity)
        real, intent(in) :: Temperature, RelativeHumidity

        T = Temperature
        RH = RelativeHumidity

        call set_coeff_drdt(T, RH)          !温湿度依存の係数の設定

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

    end function

    function read_backup(fname) result(dGroup_read)
        character(*), intent(in) :: fname
        type(dropletGroup) dGroup_read
        integer i, n_unit, num_drop
    
        print*, 'READ:', trim(fname)
        open(newunit=n_unit, form='unformatted', file=fname, status='old')
            read(n_unit) num_drop

            allocate(dGroup_read%droplet(num_drop))

            do i = 1, num_drop
                read(n_unit) dGroup_read%droplet(i)
            end do
        close(n_unit)
    
    end function read_backup
    
    function read_droplet_VTK(fname) result(dGroup_read)
        implicit none
        character(*), intent(in) :: fname
        type(dropletGroup) dGroup_read
        double precision, allocatable :: diameter(:)
        integer vn, n_unit, num_drop
        character(10) str
    
        print*, 'READ:', fname
        open(newunit=n_unit, file=fname, status='old')
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

end module dropletGroup_m