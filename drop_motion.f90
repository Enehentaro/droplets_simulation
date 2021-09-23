module drop_motion_mod
    use flow_field
    implicit none

    type :: virus_droplet
        double precision :: coordinate(3), velocity(3)
        double precision radius, radius_min
        integer :: status=0, cell_ref=0, bound_adhes=0
    end type virus_droplet
    type(virus_droplet), allocatable, private :: droplets(:)
    type(virus_droplet), allocatable, private :: droplets_ini(:)

    integer, private :: num_droplets   !全飛沫数
    integer interval, T, RH
 
    double precision dt, U_chara, L_chara, Roh_chara, Mu_chara

    double precision coeff  !蒸発方程式の係数
    double precision gumma  !空気と飛沫（水）の密度比
 
    character path_out_base*20, path_out*99, head_out*10

    double precision :: center_posi(3), width_posi(3)   !初期配置帯の中心座標および幅
    double precision :: G(3), direction_g(3)            !無次元重力加速度ベクトル, 重力方向ベクトル

    integer num_restart, n_start, n_end
    integer, private :: LoopS, LoopF, OFFSET
    integer :: num_NCS=0    !NearestCell探索を行った回数のカウンター
    double precision Rdt    !飛沫計算と気流計算の時間間隔の比
    double precision, private :: Re     !レイノルズ数

    contains

    subroutine input_condition
        double precision DTa
        integer L, n_unit

        OPEN(newunit=n_unit,FILE='condition_virus.txt',STATUS='OLD')
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
            read(n_unit,*) (center_posi(L), L=1,3)
            read(n_unit,*) (width_posi(L), L=1,3)
            read(n_unit,'()')
            read(n_unit,*) (direction_g(L), L=1,3)
            
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
            read(n_unit,*) L_chara
            read(n_unit,*) U_chara
            read(n_unit,*) Roh_chara
            read(n_unit,*) Mu_chara

        CLOSE(n_unit)

        Rdt = dt/DTa                       !データ読み込み時に時間軸合わせるパラメータ
        Re = (Roh_chara*U_chara*L_chara)/Mu_chara

        allocate(droplets(num_droplets))
        if(num_restart <= 0) allocate(droplets_ini(num_droplets))

        print*, 'restart =',num_restart
        print*, 'n_end =',n_end
        print*, 'interval =',interval
        print*, 'interval_air =',INTERVAL_FLOW
        print*, 'loop=',loops,loopf
        print*, 'dt =',dt
        print*, 'Rdt', Rdt
        print*, 'Re =',Re

        print*, 'PATH_AIR=', PATH_AIR


    end subroutine input_condition

    subroutine calc_initial_droplet
        implicit none

        if(num_restart==0) then
            call calc_initial_position
            call calc_initial_radius

        else if(num_restart==-1) then
            droplets_ini(:) = read_droplet_VTK('initial_distribution.vtk')

        else
            return  !リスタートなら無視

        end if

        droplets_ini(:)%radius_min = get_minimum_radius(droplets_ini(:)%radius) !最小半径の計算
    
    end subroutine calc_initial_droplet

    subroutine calc_initial_radius
        use csv_reader
        integer vn, i
        double precision, allocatable :: threshold(:,:)
        integer, allocatable :: rad_cnt(:)
        double precision :: radius_dim(num_droplets)
        real(8) random_rad
        !====================================================================================

        call csv_reader_dble('radius_distribution.csv', threshold)

        allocate(rad_cnt(size(threshold, dim=2)), source=0)

        call randomset  !実行時刻に応じた乱数シード設定

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

        droplets_ini(:)%radius = radius_dim(:) / L_chara         !初期飛沫半径のセットおよび無次元化

        if (sum(rad_cnt) /= num_droplets) then
            print*, 'random_rad_ERROR', sum(rad_cnt), num_droplets
            stop
        end if

        do vn = 1, size(rad_cnt)
            print*,'rad_cnt(',vn,') =',rad_cnt(vn)
        end do

    end subroutine calc_initial_radius

    function get_minimum_radius(initial_radius) result(minimum_radius)
        implicit none
        double precision, intent(in) :: initial_radius(:)
        double precision :: minimum_radius(size(initial_radius))

        if(RH==100)then
            minimum_radius(:) = initial_radius(:)
        else if(RH > 90.and.RH < 100)then
            minimum_radius(:) = initial_radius(:)*(0.0001*exp(0.0869*dble(RH)))
        else if(RH==90)then
            minimum_radius(:) = initial_radius(:)*0.28d0      !飛沫の最小半径
        else if(RH >= 64.and.RH < 90)then
            minimum_radius(:) = initial_radius(:)*(0.073*exp(0.014*dble(RH)))
        else if(RH < 64)then
            minimum_radius(:) = initial_radius(:)*0.19d0
        else
            print*,'rad_minERROR', RH
            STOP
        end if

    end function get_minimum_radius

    subroutine calc_initial_position
        integer kx,ky,kz, num_node, m, k
        double precision :: standard(3), delta(3), randble(3)

        m = 1
        k = 0

        print*, 'calc_initial_position', m, k

        standard(:) = center_posi(:) - 0.5d0*width_posi(3)

        do while(m**3 < num_droplets)
            m = m + 1
        end do
        num_node = m - 1

        delta(:) = width_posi(:) / dble(num_node-1)

        do kx = 1, num_node

            do ky = 1, num_node

                do kz = 1, num_node

                    k = k + 1
                    droplets_ini(k)%coordinate(1) = standard(1) + delta(1)*dble(kx - 1)
                    droplets_ini(k)%coordinate(2) = standard(2) + delta(2)*dble(ky - 1)
                    droplets_ini(k)%coordinate(3) = standard(3) + delta(3)*dble(kz - 1)
                    
                end do
            end do

        end do

        do while(k < num_droplets)
            k = k + 1
            call random_number(randble(:))
            droplets_ini(k)%coordinate(:) = standard(:) + width_posi(:)*randble(:)
        end do

    end subroutine calc_initial_position

    subroutine initialization_droplet
        implicit none
        character(99) fname

        if(num_restart > 0) then

            print*, 'RESTRAT'
            
            write(fname,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') num_restart
            droplets(:) = read_droplet_VTK(fname)

            block
                character(99) fname_first
                type(virus_droplet) ::  droplets_first(num_droplets)

                write(fname_first,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') 0
                droplets_first(:) = read_droplet_VTK(fname_first)
                droplets(:)%radius_min = get_minimum_radius(droplets_first(:)%radius) !最小半径の計算
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
        type(virus_droplet) :: droplets_read(num_droplets)
        double precision :: diameter(num_droplets)
        integer vn, n_unit
        !=======================================================================
    
        print*, 'READ:', fname
        open(newunit=n_unit, file=fname, status='old')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'(3(F20.16,2X))') droplets_read(vn)%coordinate(:)
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'()')
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'()')
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'(F20.16)') diameter(vn)
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'(I12)') droplets_read(vn)%status
            END DO
            read(n_unit,'()')
            DO vn = 1,num_droplets
                read(n_unit,'(3(F20.16,2X))') droplets_read(vn)%velocity(:)
            END DO
        close(n_unit)
    
        droplets_read(:)%radius = diameter(:) / 2.0d0
      
    end function read_droplet_VTK

    !***********************************************************************
    subroutine output_droplet_VTK(droplets_out, step)
        implicit none
        integer, intent(in) :: step
        type(virus_droplet), intent(in) :: droplets_out(:)
        integer vn, n_unit
        character(99) fname
        !======================================================================= 

        write(fname,'("'//trim(path_out)//trim(head_out)//'",i8.8,".vtk")') step

        !=======ここから飛沫データ（VTKファイル）の出力===========================
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ',num_droplets,' float'                              !節点の数
            DO vn = 1,num_droplets                                                             !節点の数だけループ
                write(n_unit,'(3(F20.16,2X))') droplets_out(vn)%coordinate(:)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12,2X,I12)') 'CELLS ', num_droplets, num_droplets*2                          !セルの数、セルの数×2
            DO vn = 1,num_droplets                                                             !セルの数だけループ
                write(n_unit,'(2(I12,2X))')  1, vn-1                                          !まずセルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_TYPES',num_droplets
            DO vn = 1,num_droplets                                                             !セルの数だけループ
                write(n_unit,'(I12)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_DATA ',num_droplets                                      !ここからセルのデータという合図、セルの数
            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !まずは飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1,num_droplets                                                             !セルの数だけループ
                write(n_unit,'(F20.16)') droplets_out(vn)%radius*2.0d0                            !飛沫の直径
            END DO
            write(n_unit,'(A)') 'SCALARS Status int'                                   !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1,num_droplets                                                             !セルの数だけループ
                write(n_unit,'(I12)') droplets_out(vn)%status                                            !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO
            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1,num_droplets                                                             !セルの数だけループ
                write(n_unit,'(3(F20.16,2X))') droplets_out(vn)%velocity(:)               !飛沫の速度
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
        
            write(n_unit,*) step*dt*L_chara/U_chara, ',', count(droplets_out(:)%status==0), ',',&
                count(droplets_out(:)%status==2), ',', count(droplets_out(:)%status==3)
        close(n_unit)

        ! if(count(adhesion==0) <= 0) then !浮遊粒子がなくなれば計算終了
        !   print*,'all viruses terminated',N
        !   STOP
        ! end if

    end subroutine output_droplet_VTK

    !*************************************************************************************
    subroutine output(step)
        integer, intent(in) :: step
        call output_droplet_VTK(droplets, step)
    end subroutine output

    !*************************************************************************************
    subroutine Set_Coefficients
        !*************************************************************************************
        !=====================================================================================
        double precision Es, TK, norm
        double precision, parameter :: Rv = 461.51d0                           ! 水蒸気の気体定数[J/(kg.K)]
        double precision, parameter :: Roh_d = 0.99822d0*1.0d3          ! 粒子（水）の密度[kg/m3]
        double precision, parameter :: T0 = 273.15d0                               ! [K]
        double precision, parameter :: D = 0.2564d0*1.0d-4           ! 水蒸気の拡散定数[m2/s]
        double precision, parameter :: Lv = 2.451d0*1.0d6                  ! 水の蒸発潜熱[J/kg]
        double precision, parameter :: Es0 = 6.11d0*1.0d2                  ! 基準温度における飽和蒸気圧[Pa]
        double precision, parameter :: G_dim = 9.806650d0                          ! 重力加速度[m/s2]
        !=====================================================================================  

        TK = dble(T) + T0                                    ! 室温を絶対温度[K]に変換
        Es = Es0*exp((Lv/Rv)*(1.0d0/T0 - 1.0d0/TK))       ! 室温に置ける飽和蒸気圧

        norm = norm2(direction_g(:))
        G(:) = G_dim * L_chara/(U_chara*U_chara) / norm * direction_g(:)    !無次元重力加速度
        print*, 'G_no_dim=', G(:)

        coeff = -D/(U_chara*L_chara) * (1.0d0 - dble(RH)/100.d0)*Es / (Roh_d*Rv*TK) ! dr/dt の無次元係数
        print*, 'coeff=', coeff
    
        gumma = Roh_chara / Roh_d     !密度比:    空気密度 / 飛沫(水)密度
        print*, 'gumma=', gumma

    end subroutine Set_Coefficients
    !*******************************************************************************************

    !*******************************************************************************************
    subroutine evaporation(vn) !CALCULATE drplet evaporation
        integer, intent(in) :: vn
        real(8) drdt1,R1,R_approxi,drdt2,R2
      
            !========= 飛沫半径の変化の計算　(2次精度ルンゲクッタ（ホイン）) ===========================
      
        if (droplets(vn)%radius <= droplets(vn)%radius_min) then  !半径が最小になったものを除く
            droplets(vn)%radius = droplets(vn)%radius_min
            return
        end if
    
        drdt1 = coeff / droplets(vn)%radius
        R1 = dt*drdt1
        R_approxi = droplets(vn)%radius + R1
        if(R_approxi <= 0.0d0) R_approxi = droplets(vn)%radius
        drdt2 = coeff / R_approxi
        R2 = dt*drdt2
        droplets(vn)%radius = droplets(vn)%radius + (R1+R2)*0.5d0
        
        droplets(vn)%radius = max(droplets(vn)%radius, droplets(vn)%radius_min)
      
        !*******************************************************************************************
    end subroutine evaporation
    !*******************************************************************************************

    !=====================================================================
    subroutine survival_check(step)
        integer,intent(in) :: step
        integer vfloat, vn
        double precision, save :: death_rate = 0.d0
            
        vfloat = count(droplets(:)%status == 0)
        if(vfloat == 0) return
        death_rate = death_rate + dble(vfloat)*(survival_rate(step) - survival_rate(step+1))    !このステップで死滅すべき飛沫数
        
        vn = num_droplets  !チェックはIDの後ろから
        do while(death_rate >= 1.0d0)
            if (vn==0) then
                exit

            else if (droplets(vn)%status == 0) then !浮遊粒子からのみ除去する
                droplets(vn)%status = -1
                droplets(vn)%coordinate(:) = -1.0d0
                droplets(vn)%velocity(:) = 0.0d0
                death_rate = death_rate - 1.0d0

            end if
            vn = vn - 1
        end do
      
    end subroutine survival_check

    double precision function survival_rate(step)
        integer, intent(in) :: step
        !-------- Calculate survival rate of virus ------------------------------
        if(RH == 80)then  !　相対湿度80%の時使用
            survival_rate = 0.67d0*0.5102d0**(((L_chara/U_chara)*dt*dble(step-1))/3600.0d0)
        else if(RH == 50)then  !　相対湿度50%の時使用
            survival_rate = 0.84d0*0.5735d0**(((L_chara/U_chara)*dt*dble(step-1))/3600.0d0)
        else if(RH == 35)then  !　相対湿度35%の時使用
            survival_rate = 0.86d0*0.9240d0**(((L_chara/U_chara)*dt*dble(step-1))/3600.0d0)
        else !新型コロナウイルス（1.1時間で半減）(論文によると、湿度30,60,90%のときのデータしかない)
            survival_rate = 0.999825d0**((L_chara/U_chara)*dt*dble(step-1))
        end if
    end function survival_rate

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

    subroutine motion_calc(vn)
        integer, intent(in) :: vn
        double precision  :: X(3), V(3)
        integer NCN
        double precision randble
        logical stopflag
    
        X(:) = droplets(vn)%coordinate(:)
        V(:) = droplets(vn)%velocity(:)    
        NCN = droplets(vn)%cell_ref     !前回参照セルを代入
    
        if(NCN == 0) then   !参照セルが見つかっていない（＝初期ステップ）
                NCN = nearest_cell(X)    
                print*, 'FirstNCN:', NCN
                droplets(:)%cell_ref = NCN !全粒子が同一セル参照と仮定して時間短縮を図る
    
        else
                NCN = nearer_cell(X, NCN)
                if (NCN == 0) then
                    print*, 'NCN_ERROR:', vn, NCN
                    stop
                end if
                if (.not.nearcell_check(X(:), NCN)) NCN = nearest_cell(X)
    
        end if

        stopflag = .false.
        if (NoB(NCN) >= 1) stopflag = adhesion_check(vn, NCN)
    
        droplets(vn)%cell_ref = NCN    !結果を参照セル配列に記憶

        call area_check(vn, stopflag)
    
        if (stopflag) then
            
                ! if ((X(1)>1.59d0.and.X(1)<2.01d0).and.&
                !       (X(2)>0.25d0.and.X(2)<0.35d0).and.&
                !       (X(3)>0.01d0.and.X(3)<0.75d0)) then
                !       droplets(vn)%status = 2       !AP_air_cleaner
                ! else
                !       droplets(vn)%status = 1
                ! end if
    
                if ((X(1)>0.29d0.and.X(1)<0.31d0).and.&
                        (X(2)>2.94d0.and.X(2)<3.36d0).and.&
                        (X(3)>0.01d0.and.X(3)<0.74d0)) then
    
                    droplets(vn)%status = 2       !ACAP_air_cleaner_left
    
                else if ((X(1)>1.58d0.and.X(1)<2.02d0).and.&
                        (X(2)>5.99d0.and.X(2)<6.01d0).and.&
                        (X(3)>0.01d0.and.X(3)<0.74d0)) then
    
                    droplets(vn)%status = 2       !ACAP_air_cleaner_oposit
    
                else if ((X(1)>1.58d0.and.X(1)<2.02d0).and.&
                        (X(2)>0.29d0.and.X(2)<0.31d0).and.&
                        (X(3)>0.01d0.and.X(3)<0.74d0)) then
    
                    droplets(vn)%status = 2       !ACAP_air_cleaner_under
    
    
                else if ((X(1)>1.43d0.and.X(1)<2.18d0).and.&
                        (X(2)>0.01d0.and.X(2)<0.22d0).and.&
                        (X(3)>2.15d0.and.X(3)<2.18d0)) then
    
                    call random_number(randble)
    
                    if (randble < 0.015d0) then
                            droplets(vn)%status = 3       !ACAP_air_conditioner
    
                    else
                            droplets(vn)%coordinate(1) = X(1)
                            droplets(vn)%coordinate(2) = X(2) + 0.23d0
                            droplets(vn)%coordinate(3) = X(3) - 0.26d0
    
                            droplets(vn)%velocity(:) = 0.0d0     !速度をゼロに
    
                            return
            
                    end if
    
    
                else
                    droplets(vn)%status = 1

                end if
    
                droplets(vn)%velocity(:) = 0.0d0     !速度をゼロに
                droplets(vn)%coordinate(:) = X(:)
    
        else
    
                droplets(vn)%velocity(:) = get_velocity(V(:), VELC(:, NCN), droplets(vn)%radius)
                
                droplets(vn)%coordinate(:) = X(:) + (V(:) + droplets(vn)%velocity(:))*0.5d0*dt
            
        end if
    
        
    end subroutine motion_calc

    !*******************************************************************************************
    function get_velocity(vel_d, vel_a, radius_d) result(vel_d_next)
        !*******************************************************************************************
        !=====================================================================================
        double precision, intent(in) :: vel_d(3), vel_a(3), radius_d
        double precision speed_r, Re_d, Cd, Coefficient, vel_d_next(3)
        !=====================================================================================

        speed_r = norm2(vel_a(:) - vel_d(:))
        Re_d = (speed_r * 2.0d0*radius_d) * Re + 1.d-9  !ゼロ割回避のため、小さな値を足す

        Cd = (24.0d0/Re_d)*(1.0d0 + 0.15d0*(Re_d**0.687d0))

        Coefficient = (3.0d0*Cd*gumma*speed_r)/(8.0d0*radius_d)

        vel_d_next(:) = ( vel_d(:) + ( G(:) + Coefficient*vel_a(:) )*dt ) &
                            / ( 1.0d0 + Coefficient*dt )

    end function get_velocity
    !----------------------------------------------------------------------------------

!*******************************************************************************************
    integer function nearest_cell(X) !粒子vnに最も近いセルNCNの探索
        double precision, intent(in) :: X(3)
        integer II, IIMX
        double precision, allocatable :: distance(:)
        !=====================================================================================
        num_NCS = num_NCS +1

        IIMX = size(CENC, dim=2)

        allocate(distance(IIMX))
        !↓↓↓↓　一番近いセル中心の探索
        !$omp parallel do
        DO II = 1,IIMX
                distance(II) = norm2(CENC(:,II) - X(:))
        END DO
        !$omp end parallel do 
        !↑↑↑↑
        
        nearest_cell = minloc(distance, dim=1)   !最小値インデックス
        
    end function nearest_cell

!----------------------------------------------------------------------------------   
!*******************************************************************************************
    integer function nearer_cell(X, NCN)  !近セルの探索（隣接セルから）
        integer, intent(in) :: NCN
        double precision, intent(in) :: X(3)
        integer NC, IIaround, index_min
        double precision :: distancecheck(2)
        double precision, allocatable :: distance(:)
        !=====================================================================================
        nearer_cell = NCN
        allocate(distance(NCMAX))
        distancecheck(1) = norm2(CENC(:,nearer_cell)-X(:))   !注目セル重心と粒子との距離
        
        check:DO
                distance(:) = 1.0d10     !初期値はなるべく大きくとる
        
                DO NC = 1, NUM_NC(nearer_cell)  !全隣接セルに対してループ。
                    IIaround = NEXT_CELL(NC, nearer_cell)       !現時点で近いとされるセルの隣接セルのひとつに注目
                    IF (IIaround > 0) then
                            distance(NC) = norm2(CENC(:,IIaround)-X(:))   !注目セル重心と粒子との距離を距離配列に代入
                    END IF
        
                END DO
        
                distancecheck(2) = minval(distance,dim=1)     !距離配列の最小値
        
                if(distancecheck(2) < distancecheck(1)) then !より近いセルの発見で条件満足
                    distancecheck(1) = distancecheck(2)    !最小値の更新
                    index_min = minloc(distance,dim=1)            !最小値のインデックス
                    nearer_cell = NEXT_CELL(index_min, nearer_cell)    !現時点で近いとされるセルの更新
                    if(nearer_cell==0) then
                            print*,'nearer_cell_error', nearer_cell, X(:)
                            return
                    end if

                else  !より近いセルを発見できなかった場合
        
                    exit check     !ループ脱出
        
                end if
        
        END DO check
        
    end function nearer_cell
!*******************************************************************************************
                    
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

                    
    logical function nearcell_check(X, NCN)
        double precision, intent(in) :: X(3)
        integer, intent(in) :: NCN
        double precision :: distance

        distance = norm2(X(:)-CENC(:,NCN))

        if (distance < 1.0d1*WIDC(NCN)) then
                nearcell_check = .True.
        else
                nearcell_check = .False.
        end if


    end function nearcell_check


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

        write(digits_fmt,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS

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

        do vn = 1, size(droplets)
        
                if (droplets(vn)%status <= 0) cycle !付着していないならスルー
        
                JB = droplets(vn)%bound_adhes
                if (JB > 0) then
                    droplets(vn)%coordinate(:) = droplets(vn)%coordinate(:) + CENF(:,JB,2) - CENF(:,JB,1) !面重心の移動量と同じだけ移動
                else
                    call area_check(vn)
        
                end if
        
        end do

        ! print*, 'FIN:boundary_move'

    end subroutine boundary_move

    subroutine area_check(vn,check)
        integer, intent(in) :: vn
        logical,optional,intent(inout) :: check
        integer L

        do L = 1, 3
    
            if(droplets(vn)%coordinate(L) < MIN_CDN(L)) then
                droplets(vn)%coordinate(L) = MIN_CDN(L)
                if(present(check)) check = .true.
            else if(droplets(vn)%coordinate(L) > MAX_CDN(L)) then
                droplets(vn)%coordinate(L) = MAX_CDN(L)
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

            case default
                get_drop_info = 0

        end select

    end function get_drop_info

    subroutine randomset    !実行時刻に依存した乱数シードを指定する
        implicit none
        integer :: seedsize, i
        integer, allocatable :: seed(:)

        print*, 'call:randomset'
    
        call random_seed(size=seedsize) !シードのサイズを取得。（コンパイラごとに異なるらしい）
        allocate(seed(seedsize)) !新シード配列サイズの割り当て
    
        do i = 1, seedsize
            call system_clock(count=seed(i)) !時間を新シード配列に取得
        end do
    
        call random_seed(put=seed(:)) !新シードを指定
          
    end subroutine randomset


end module drop_motion_mod