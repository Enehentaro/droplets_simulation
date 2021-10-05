module drop_motion_mod
    use flow_field
    use equation_mod
    implicit none

    type :: virus_droplet
        double precision :: coordinate(3), velocity(3)=0.d0
        double precision radius, radius_min
        integer :: status=0, cell_ref=0, bound_adhes=0
    end type virus_droplet
    type(virus_droplet), allocatable, private :: droplets(:)
    type(virus_droplet), allocatable, private :: droplets_ini(:)

    type, extends(virus_droplet) :: vd_extracted
        integer original_ID
    end type

    integer, private :: num_droplets   !全飛沫数
    integer interval, T, RH
 
    character path_out_base*99, path_out*99, head_out*10

    integer num_restart, n_start, n_end
    integer, private :: LoopS, LoopF, OFFSET
    double precision Rdt    !飛沫計算と気流計算の時間間隔の比

    contains

    subroutine pre_setting

        call read_and_set_condition

        if(num_restart==0) then
            
            call random_set  !実行時刻に応じた乱数シード設定
            call calc_initial_position
            call calc_initial_radius

        else if(num_restart==-1) then
            droplets_ini(:) = read_droplet_VTK('initial_distribution.vtk')

        else
            return  !リスタートなら無視

        end if

        droplets_ini(:)%radius_min = get_minimum_radius(droplets_ini(:)%radius, RH) !最小半径の計算

    end subroutine pre_setting

    subroutine read_and_set_condition
        double precision DTa, dt, L, U, Rho, Mu
        double precision :: direction_g(3)
        integer i, n_unit

        OPEN(newunit=n_unit,FILE='condition.txt',STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,*) num_restart
            read(n_unit,'()')
            read(n_unit,*) n_end
            read(n_unit,'()')
            read(n_unit,*) dt
            read(n_unit,'()')
            read(n_unit,*) interval
            read(n_unit,'()')
            read(n_unit,'(A)') path_out_base
            read(n_unit,*) head_out
            read(n_unit,'()')
            read(n_unit,*) T
            read(n_unit,*) RH
            read(n_unit,'()')
            read(n_unit,*) num_droplets
            read(n_unit,'()')
            read(n_unit,*) (direction_g(i), i=1,3)
            
            read(n_unit,'()')
    
            read(n_unit,'()')
            read(n_unit,'(A)') PATH_AIR
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
            read(n_unit,*) Rho
            read(n_unit,*) Mu

        CLOSE(n_unit)

        Rdt = dt/DTa                       !データ読み込み時に時間軸合わせるパラメータ

        if(num_restart <= 0) then   !リスタートでないなら割り付け
            allocate(droplets(num_droplets))
            allocate(droplets_ini(num_droplets))
        end if

        print*, 'restart =',num_restart
        print*, 'n_end =',n_end
        print*, 'interval =',interval
        print*, 'interval_air =',INTERVAL_FLOW
        print*, 'loop=',loops,loopf
        print*, 'dt =',dt
        print*, 'Rdt', Rdt

        print*, 'PATH_AIR=', PATH_AIR

        call set_basical_variables(dt, L, U, Rho, Mu)

        call set_gravity_acceleration(direction_g)

    end subroutine read_and_set_condition

    subroutine calc_initial_radius
        use csv_reader
        integer vn, i
        double precision, allocatable :: threshold(:,:)
        integer, allocatable :: rad_cnt(:)
        double precision :: radius_dim(num_droplets)
        real(8) random_rad
        !====================================================================================

        ! call csv_reader_dble('radius_distribution.csv', threshold)
        threshold = read_csv_dble('radius_distribution.csv')

        allocate(rad_cnt(size(threshold, dim=2)), source=0)

        do vn = 1, num_droplets                       !飛沫半径の分布を乱数によって与える

            call random_number(random_rad)

            do i = 1, size(threshold, dim=2)
                if(random_rad < threshold(2, i)) then
                    radius_dim(vn) = threshold(1, i) * 1.0d-6
                    rad_cnt(i) = rad_cnt(i) + 1
                    exit
                end if
            end do

        end do

        print*, shape(droplets_ini(:)%radius), shape(radius_dim(:))

        droplets_ini(:)%radius = radius_dim(:) / representative_value('Length') !初期飛沫半径のセットおよび無次元化

        if (sum(rad_cnt) /= num_droplets) then
            print*, 'random_rad_ERROR', sum(rad_cnt), num_droplets
            stop
        end if

        do vn = 1, size(rad_cnt)
            print*,'rad_cnt(',vn,') =',rad_cnt(vn)
        end do

    end subroutine calc_initial_radius

    subroutine calc_initial_position
        use csv_reader
        integer kx,ky,kz, num_per_edge, num_per_box, m, k, k_end, cnt
        integer i_box, num_box
        double precision :: standard(3), delta(3), width(3), randble(3)
        double precision, allocatable :: position_mat(:,:)

        call csv_reader_dble('initial_position.csv', position_mat)

        num_box = size(position_mat, dim=2)

        num_per_box = num_droplets / num_box
        
        print*, 'calc_initial_position'

        m = 1
        do while(num_box*(m**3) < num_droplets)
            m = m + 1
        end do
        num_per_edge = m - 1    !配置帯一辺当たりの飛沫数

        k = 0
        cnt = 0
        do i_box = 1, num_box
  
            width(:) = position_mat(4:6, i_box)
            standard(:) = position_mat(1:3, i_box) - 0.5d0*width(:)
            delta(:) = width(:) / dble(num_per_edge - 1)

            do kx = 1, num_per_edge

                do ky = 1, num_per_edge

                    do kz = 1, num_per_edge

                        k = k + 1
                        droplets_ini(k)%coordinate(1) = standard(1) + delta(1)*dble(kx - 1)
                        droplets_ini(k)%coordinate(2) = standard(2) + delta(2)*dble(ky - 1)
                        droplets_ini(k)%coordinate(3) = standard(3) + delta(3)*dble(kz - 1)
                        
                    end do

                end do

            end do

            k_end = i_box * num_per_box
            if(i_box == num_box) k_end = num_droplets

            do while(k < k_end)
                k = k + 1
                call random_number(randble(:))
                droplets_ini(k)%coordinate(:) = standard(:) + width(:)*randble(:)
            end do

            print*, 'BOX', i_box, 'has', k - cnt, 'droplets.'

            cnt = k

        end do

    end subroutine calc_initial_position

    subroutine initialization_droplet
        implicit none
        character(99) fname

        if(num_restart > 0) then

            print*, 'RESTRAT'
            
            write(fname,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') num_restart
            droplets = read_droplet_VTK(fname)   !ここで自動割り付け
            num_droplets = size(droplets(:))

            block
                character(99) fname_first
                type(virus_droplet), allocatable ::  droplets_first(:)

                write(fname_first,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') 0
                droplets_first = read_droplet_VTK(fname_first)   !自動割り付け
                droplets(:)%radius_min = get_minimum_radius(droplets_first(:)%radius, RH) !最小半径の計算
            end block

            n_start = num_restart

        else
            droplets(:) = droplets_ini(:)
            n_start = 0
            call output(n_start)  !リスタートでないなら初期配置出力

        end if
            
    end subroutine initialization_droplet

    !***********************************************************************
    function read_droplet_VTK(fname) result(droplets_read)
        implicit none
        character(*), intent(in) :: fname
        type(virus_droplet), allocatable :: droplets_read(:)
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

            allocate(droplets_read(num_drop), diameter(num_drop))

            DO vn = 1, num_drop
                read(n_unit, *) droplets_read(vn)%coordinate(:)
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
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) diameter(vn)
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit,'(I12)') droplets_read(vn)%status
            END DO
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets_read(vn)%velocity(:)
            END DO
        close(n_unit)
    
        droplets_read(:)%radius = diameter(:) * 0.5d0
      
    end function read_droplet_VTK

    !***********************************************************************
    subroutine output_droplet_VTK(droplets_out, step)
        implicit none
        integer, intent(in) :: step
        type(virus_droplet), intent(in) :: droplets_out(:)
        integer vn, n_unit, num_drop
        character(99) fname
        !======================================================================= 

        num_drop = size(droplets_out(:))
        write(fname,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') step

        !=======ここから飛沫データ（VTKファイル）の出力===========================
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ',num_droplets,' float'                              !節点の数
            DO vn = 1, num_drop                                                            !節点の数だけループ
                write(n_unit,'(3(E20.12e2,2X))') droplets_out(vn)%coordinate(:)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12,2X,I12)') 'CELLS ', num_droplets, num_droplets*2                          !セルの数、セルの数×2
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(2(I12,2X))')  1, vn-1                                          !まずセルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_TYPES',num_droplets
            DO vn = 1, num_drop                                                           !セルの数だけループ
                write(n_unit,'(I12)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_DATA ',num_droplets                                      !ここからセルのデータという合図、セルの数
            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !まずは飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(E20.12e2)') droplets_out(vn)%radius*2.0d0                            !飛沫の直径
            END DO
            write(n_unit,'(A)') 'SCALARS Status int'                                   !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I12)') droplets_out(vn)%status                                            !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO
            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1, num_drop                                                             !セルの数だけループ
                write(n_unit,'(3(E20.12e2,2X))') droplets_out(vn)%velocity(:)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'WRITEOUT:', fname

    !=======飛沫データ（VTKファイル）の出力ここまで===========================

        !以下はCSVファイルの出力

        if(step==0) then !初期ステップならファイル新規作成
            open(newunit=n_unit, file=trim(path_out)//'particle_'//trim(head_out)//'.csv', status='replace')
            print*,'REPLACE:particle_data.csv'

        else
            open(newunit=n_unit, file=trim(path_out)//'particle_'//trim(head_out)//'.csv'&
                , action='write', status='old', position='append')

        end if
        
            write(n_unit,*) real_time(step), ',', count(droplets_out(:)%status==0), ',',&
                count(droplets_out(:)%status==2), ',', count(droplets_out(:)%status==3)
        close(n_unit)

        ! if(count(adhesion==0) <= 0) then !浮遊粒子がなくなれば計算終了
        !   print*,'all viruses terminated',N
        !   STOP
        ! end if

    end subroutine output_droplet_VTK

    subroutine output_droplet_CSV(droplets_out, step)
        implicit none
        integer, intent(in) :: step
        type(virus_droplet), intent(in) :: droplets_out(:)
        integer n_unit

        !以下はCSVファイルの出力

        if(step==0) then !初期ステップならファイル新規作成
            open(newunit=n_unit, file=trim(path_out)//'particle_'//trim(head_out)//'.csv', status='replace')
            print*,'REPLACE:particle_data.csv'

        else
            open(newunit=n_unit, file=trim(path_out)//'particle_'//trim(head_out)//'.csv'&
                , action='write', status='old', position='append')

        end if
        
            write(n_unit,*) real_time(step), ',', count(droplets_out(:)%status==0), ',',&
                count(droplets_out(:)%status==2), ',', count(droplets_out(:)%status==3)
        close(n_unit)

        ! if(count(adhesion==0) <= 0) then !浮遊粒子がなくなれば計算終了
        !   print*,'all viruses terminated',N
        !   STOP
        ! end if

    end subroutine output_droplet_CSV

    subroutine output(step)
        integer, intent(in) :: step
        call output_droplet_VTK(droplets, step)
        call output_droplet_CSV(droplets, step)
    end subroutine output

    !=====================================================================
    subroutine survival_check(step)
        integer,intent(in) :: step
        integer vfloat, vn
        double precision rand
        double precision, save :: death_rate = 0.d0
            
        vfloat = count(droplets(:)%status == 0)
        if(vfloat == 0) return  !浮遊数がゼロならリターン
        death_rate = death_rate + dble(vfloat)*(survival_rate(step) - survival_rate(step+1))    !このステップで死滅すべき飛沫数
        
        do while(death_rate >= 1.0d0)
            call random_number(rand)    !死滅IDは乱数で決まる
            vn = int(num_droplets*rand)
            if (droplets(vn)%status == 0) then !浮遊粒子からのみ除去する
                droplets(vn)%status = -1
                droplets(vn)%coordinate(:) = MIN_CDN(:) - 1.0d0 !計算エリア外に配置（不要かも）
                droplets(vn)%velocity(:) = 0.0d0
                death_rate = death_rate - 1.0d0

            end if
        end do
      
    end subroutine survival_check

    subroutine Calculation_Droplets
        implicit none
        integer vn
        
        !$omp parallel do private(vn)
        DO vn = 1, num_droplets
            if(droplets(vn)%status /= 0) cycle !浮遊状態でないなら無視
            call evaporation(vn)    !半径変化方程式
            call motion_calc(vn)     !運動方程式
        END DO
        !$omp end parallel do 

    end subroutine Calculation_Droplets

        !*******************************************************************************************
    subroutine evaporation(vn) !CALCULATE drplet evaporation
        integer, intent(in) :: vn
        double precision radius_n
      
            !========= 飛沫半径の変化の計算　(2次精度ルンゲクッタ（ホイン）) ===========================
      
        if (droplets(vn)%radius <= droplets(vn)%radius_min) then  !半径が最小になったものを除く
            droplets(vn)%radius = droplets(vn)%radius_min
            return
        end if
    
        radius_n = evaporatin_eq(droplets(vn)%radius)
        
        droplets(vn)%radius = max(radius_n, droplets(vn)%radius_min)
      
        !*******************************************************************************************
    end subroutine evaporation
    !*******************************************************************************************

    subroutine motion_calc(vn)
        integer, intent(in) :: vn
        double precision  :: X(3), V(3)
        integer NCN
        logical stopflag
    
        X(:) = droplets(vn)%coordinate(:)
        V(:) = droplets(vn)%velocity(:)

        droplets(vn)%cell_ref = search_cell_ref(X(:), droplets(vn)%cell_ref)

        NCN = droplets(vn)%cell_ref

        stopflag = .false.
        if (NoB(NCN) >= 1) stopflag = adhesion_check(vn, NCN)

        call area_check(X(:), stopflag)
    
        if (stopflag) then
            droplets(vn)%status = 1
            droplets(vn)%velocity(:) = 0.0d0     !速度をゼロに
            droplets(vn)%coordinate(:) = X(:)
    
        else
    
            droplets(vn)%velocity(:) = motion_eq(V(:), VELC(:, NCN), droplets(vn)%radius)
            
            droplets(vn)%coordinate(:) = next_position(X(:), V(:), droplets(vn)%velocity(:))
            
        end if
        
    end subroutine motion_calc

    function search_cell_ref(X, NCN) result(cell_reference)
        double precision, intent(in) :: X(3)
        integer, intent(in) :: NCN
        integer cell_reference
            
        if(NCN == 0) then   !参照セルが見つかっていない（＝初期ステップ）
            cell_reference = nearest_cell(X)    
            print*, 'FirstNCN:', cell_reference
            droplets(:)%cell_ref = cell_reference !全粒子が同一セル参照と仮定して時間短縮を図る
    
        else
            cell_reference = nearer_cell(X, NCN)
            if (cell_reference == 0) then
                print*, 'NCN_ERROR:', X(:), cell_reference
                stop
            end if

            if (.not.nearcell_check(X(:), cell_reference)) cell_reference = nearest_cell(X)
    
        end if

    end function search_cell_ref
                    
    logical function adhesion_check(vn, NCN)

        integer JJ, JB
        integer, intent(in) :: vn, NCN
        double precision :: r_vector(3), inner

        adhesion_check = .false.

        do JJ = 1, NoB(NCN)
                JB = ICB(JJ, NCN)

                r_vector(:) = droplets(vn)%coordinate(:) - CENF(:,JB,1)

                inner = sum(r_vector(:)*NVECF(:,JB))

                if (inner >= 0.0d0) then
                    adhesion_check = .true. !外向き法線ベクトルと位置ベクトルの内積が正なら付着判定
                    droplets(vn)%bound_adhes = JB     !付着した境界面番号
                end if
        end do

    end function adhesion_check

    integer function get_num_air(n_virus)
        integer, intent(in) :: n_virus
        integer Lamda

        get_num_air = int(dble(N_virus)*RDT)    !気流計算における経過ステップ数に相当

        Lamda = LoopF - LoopS
        
        if((Lamda>0).and.(get_num_air>(LoopF-OFFSET))) get_num_air = mod(get_num_air, Lamda)

        get_num_air = get_num_air + OFFSET

    end function

    subroutine read_flow_field(n_virus)
        integer, intent(in) :: n_virus
        integer FNUM
        character(99) FNAME
        character(4) digits_fmt

        FNUM = get_num_air(n_virus)

        digits_fmt = get_digits_format()

        select case(FILE_TYPE)
            case('VTK')
                if (INTERVAL_FLOW == -1) then !定常解析
                    FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
                else
                    write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".vtk")') FNUM

                end if
                call read_VTK(FNAME)

            case('INP')
                if(INTERVAL_FLOW==-1) then
                    FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
                else
                    if(FNUM==0) then
                        write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".inp")') 1
                    else
                        write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".inp")') FNUM
                    end if
                end if
                call read_INP(FNAME)   !INPを読み込む(SHARP用)

            case default
                print*,'FILE_TYPE NG:', FILE_TYPE
                STOP
                    
        end select
            
        MAX_CDN(1) = maxval(CDN(1,:))
        MAX_CDN(2) = maxval(CDN(2,:))
        MAX_CDN(3) = maxval(CDN(3,:))
        print*, 'MAX_coordinates=', MAX_CDN(:)
            
        MIN_CDN(1) = minval(CDN(1,:))
        MIN_CDN(2) = minval(CDN(2,:))
        MIN_CDN(3) = minval(CDN(3,:))
        print*, 'MIN_coordinates=', MIN_CDN(:)
            
        call set_gravity_center
        call boundary_set
        if(n_virus > n_start) call boundary_move

        CENF(:,:,1) = CENF(:,:,2)
            
    end subroutine read_flow_field
                      
    subroutine boundary_move !境界面の移動に合わせて付着飛沫も移動
        integer vn, JB

        ! print*, 'CALL:boundary_move'

        do vn = 1, num_droplets
        
            if (droplets(vn)%status <= 0) cycle !付着していないならスルー
    
            JB = droplets(vn)%bound_adhes
            if (JB > 0) then
                droplets(vn)%coordinate(:) = droplets(vn)%coordinate(:) + CENF(:,JB,2) - CENF(:,JB,1) !面重心の移動量と同じだけ移動
            else
                call area_check(droplets(vn)%coordinate(:))
    
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

    integer function get_drop_info(name)
        character(*), intent(in) :: name

        select case(name)
            case('adhesion')
                get_drop_info = count(droplets(:)%status > 0)

            case('floating')
                get_drop_info = count(droplets(:)%status == 0)

            case('death')
                get_drop_info = count(droplets(:)%status == -1)

            case('coalescence')
                get_drop_info = count(droplets(:)%status == -2)

            case default
                get_drop_info = 0

        end select

    end function get_drop_info

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

    function get_floating_droplets(droplets_in) result(floating_droplets)
        implicit none
        type(virus_droplet), intent(in) :: droplets_in(:)
        type(vd_extracted), allocatable :: floating_droplets(:)
        integer i, j, num_floating

        num_floating = count(droplets_in(:)%status==0)
        allocate(floating_droplets(num_floating))

        j = 0
        do i =  1, size(droplets_in(:))
            if(droplets_in(i)%status==0) then
                j = j + 1
                floating_droplets(j)%virus_droplet = droplets_in(i)
                floating_droplets(j)%original_ID = i
            end if
        end do

    end function get_floating_droplets

    subroutine coalescence_check(step)
        integer,intent(in) :: step
        integer d1, d2, num_floating, i, i_origin
        integer, save :: last_coalescence
        type(vd_extracted), allocatable :: floatings(:)
        double precision :: distance, r1, r2

        !最後の合体から100ステップが経過したら、以降は合体が起こらないとみなしてリターン
        if((step - last_coalescence) > 100) return

        floatings = get_floating_droplets(droplets(:))
        num_floating = size(floatings(:))

        drop1 : do d1 = 1, num_floating - 1
            if(floatings(d1)%status/=0) cycle drop1

            drop2 : do d2 = d1 + 1, num_floating
                if(floatings(d2)%status/=0) cycle drop2

                distance = norm2(floatings(d2)%coordinate(:) - floatings(d1)%coordinate(:))
                r1 = floatings(d1)%radius
                r2 = floatings(d2)%radius

                if((r1+r2) >= distance) then
                    print*, 'Coalescence', d1, d2
                    call coalescence(floatings(d1)%virus_droplet, floatings(d2)%virus_droplet)
                    last_coalescence = step

                end if

            end do drop2

        end do drop1

        do i = 1, num_floating
            i_origin = floatings(i)%original_ID
            droplets(i_origin) = floatings(i)%virus_droplet
        end do

    end subroutine coalescence_check

    subroutine coalescence(droplet1, droplet2)
        type(virus_droplet), intent(inout) :: droplet1, droplet2
        double precision :: radius_c, volume1, volume2, velocity_c(3)

        volume1 = droplet1%radius**3
        volume2 = droplet2%radius**3
        radius_c = (volume1 + volume2)**(1.d0/3.d0) !結合後の飛沫半径
        velocity_c(:) = (volume1*droplet1%velocity(:) + volume2*droplet2%velocity(:)) / (volume1 + volume2)
        
        if(volume1 >= volume2) then
            droplet1%radius = radius_c
            droplet1%velocity(:) = velocity_c(:)
            droplet2%radius = 0.d0
            droplet2%velocity(:) = 0.d0
            droplet2%status = -2
        else
            droplet2%radius = radius_c
            droplet2%velocity(:) = velocity_c(:)
            droplet1%radius = 0.d0
            droplet1%velocity(:) = 0.d0
            droplet1%status = -2
        end if

    end subroutine coalescence

end module drop_motion_mod