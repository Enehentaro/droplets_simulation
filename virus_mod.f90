module virus_mod
    implicit none
    integer, private :: vnmax   !全飛沫数
    integer interval, T, RH
 
    double precision dt, U_chara, L_chara, Roh_chara, Mu_chara

    double precision coeff  !蒸発方程式の係数
    double precision gumma  !空気と飛沫（水）の密度比
 
    character path_out_base*20, path_out*99, head_out*10

    double precision :: center_posi(3), width_posi(3)   !初期配置帯の中心座標および幅
    double precision :: G(3), direction_g(3)            !無次元重力加速度ベクトル, 重力方向ベクトル
 
    integer, allocatable :: adhesion(:)                 !飛沫状態（ゼロなら浮遊、非ゼロで静止）
    integer, allocatable :: vn_trans(:)                 !浮遊飛沫ID変換
    integer, allocatable :: rad_cnt(:)                  !初期飛沫半径分布

    double precision, allocatable :: crd_drp(:,:,:)     !飛沫座標
    double precision, allocatable :: vel_drp(:,:,:)     !飛沫速度
    double precision, allocatable :: radius(:,:), rad_min(:), rad_ini(:)    !飛沫半径、最小半径、初期半径

    contains

!*************************************************************************************
    subroutine initial_virus(num_rst)
        use csv_reader
        integer,intent(in) :: num_rst
        integer vn, i, n_unit
        double precision, allocatable :: threshold(:,:)
        real(8) random_rad
        character fname*99
        !====================================================================================
        !========= Initial Position of drplets ===========================

        call set_initial_position

        vel_drp(:,:,:) = 0.0d0
        adhesion(:) = 0

        !========= drplet radius ===========================
        if(.not.allocated(rad_cnt)) then

            call csv_reader_dble('radius_distribution.csv', threshold, 3)

            allocate(rad_cnt(size(threshold, dim=2)))

            fname = 'initial_distribution.txt'

            if(num_rst==0) then !ゼロなら乱数を用いて分布決定

                rad_cnt(:) = 0

                call randomset  !実行時刻に応じた乱数シード設定

                do vn = 1, vnmax                       !飛沫半径の分布を乱数によって与える

                    call random_number(random_rad)

                    do i = 1, size(threshold, dim=2)
                        if(random_rad < threshold(2, i)) then
                            rad_ini(vn) = threshold(1, i) * 1.0d-6
                            rad_cnt(i) = rad_cnt(i) + 1
                            exit
                        end if
                    end do

                end do
                
                !========= Output initial distribution ===========================
                open(newunit=n_unit,file = trim(path_out)//'/'//fname, STATUS='REPLACE')
                    do i = 1, size(rad_cnt)
                        write(n_unit,*) rad_cnt(i)      !ウイルス半径ごとの個数
                    end do
                    do vn = 1, vnmax
                        write(n_unit,'(F20.16)') rad_ini(vn)
                    end do
                close(n_unit)
                !=================================================================

            else    !テキストファイルから分布指定

                if(num_rst >= 1) fname = trim(path_out)//'/'//fname
                !========= input initial distribution ===========================
                print*, 'READ:', fname
                open(newunit=n_unit, file=fname, STATUS='OLD')
                    do i = 1, size(rad_cnt)
                        read(n_unit,*) rad_cnt(i)      !ウイルス半径ごとの個数
                    end do
                    do vn = 1, vnmax
                        read(n_unit,*) rad_ini(vn)
                    end do
                close(n_unit)
                !=================================================================

            end if

            if(num_rst <= 0) then

                radius(:,1) = rad_ini(:) / L_chara         !初期飛沫半径のセットおよび無次元化

            else !1以上ならリスタート

                call restartREAD(num_rst)

            end if

            if (sum(rad_cnt) /= vnmax) then
                print*, 'random_rad_ERROR', sum(rad_cnt), vnmax
                stop
            end if

        end if

        do vn = 1, size(rad_cnt)
            print*,'rad_cnt(',vn,') =',rad_cnt(vn)
        end do

        !========= Set minimum radius ===========================
        if(RH==100)then
            rad_min(:) = rad_ini(:)
        else if(RH > 90.and.RH < 100)then
            rad_min(:) = rad_ini(:)*(0.0001*exp(0.0869*dble(RH)))
        else if(RH==90)then
            rad_min(:) = rad_ini(:)*0.28d0      !飛沫の最小半径
        else if(RH >= 64.and.RH < 90)then
            rad_min(:) = rad_ini(:)*(0.073*exp(0.014*dble(RH)))
        else if(RH < 64)then
            rad_min(:) = rad_ini(:)*0.19d0
        else
            print*,'rad_minERROR',RH
            STOP
        end if
        rad_min(:) = rad_min(:) / L_chara     !最小半径の無次元化
        !=================================================================

    end subroutine initial_virus


    subroutine set_initial_position

        integer kx,ky,kz, num_node
        integer ::  m=1, vn=0
        double precision :: standard(3), delta(3), randble(3)

        standard(:) = center_posi(:) - 0.5d0*width_posi(3)

        do while(m**3 < vnmax)
            m = m + 1
        end do
        num_node = m - 1

        delta(:) = width_posi(:) / dble(num_node-1)

            do kx = 1, num_node

                do ky = 1, num_node

                    do kz = 1, num_node

                        vn = vn + 1
                        crd_drp(1,vn,1) = standard(1) + delta(1)*dble(kx - 1)
                        crd_drp(2,vn,1) = standard(2) + delta(2)*dble(ky - 1)
                        crd_drp(3,vn,1) = standard(3) + delta(3)*dble(kz - 1)
                        
                    end do
                end do

            end do

            do while(vn < vnmax)
                vn = vn + 1
                call random_number(randble(:))
                crd_drp(:,vn,1) = standard(:) + width_posi(:)*randble(:)
            end do

    end subroutine set_initial_position

    !***********************************************************************
    subroutine restartREAD(step)
        integer, intent(in) :: step
        integer vn, n_unit
        character(len=99) :: FN
        !=======================================================================

        write(FN,'("'//trim(path_out)//'/'//trim(head_out)//'",i8.8,".vtk")') step
    
        print*, 'restartREAD:', FN,', step=',step
        open(newunit=n_unit, FILE=FN,STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'(3(F20.16,2X))')  crd_drp(1,vn,1),crd_drp(2,vn,1),crd_drp(3,vn,1) 
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'()')
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'()')
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'(F20.16)') radius(vn,1)
            END DO
            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'(I12)') adhesion(vn)
            END DO
            read(n_unit,'()')
            DO vn = 1,vnmax
                read(n_unit,'(3(F20.16,2X))') vel_drp(1,vn,1),vel_drp(2,vn,1),vel_drp(3,vn,1)
            END DO
        close(n_unit)
    
        radius(:,1) = radius(:,1) / 2.0d0

        crd_drp(:,:,2) = crd_drp(:,:,1)
        radius(:,2) = radius(:,1)
        vel_drp(:,:,2) = vel_drp(:,:,1)
      
        end subroutine restartREAD


    !***********************************************************************
    subroutine writeout(step)
        integer, intent(in) :: step
        integer vn, vntotal, n_unit
        character(len=50) :: FN
        !======================================================================= 

        write(FN,'("'//trim(path_out)//'/'//trim(head_out)//'",i8.8,".vtk")') step

        !=======ここから飛沫データ（VTKファイル）の出力===========================

        vntotal = vnmax*2  
        open(newunit=n_unit, FILE=FN, STATUS='REPLACE')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ',vnmax,' float'                              !節点の数
            DO vn = 1,vnmax                                                             !節点の数だけループ
                write(n_unit,'(3(F20.16,2X))') crd_drp(1,vn,1),crd_drp(2,vn,1),crd_drp(3,vn,1)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12,2X,I12)') 'CELLS ',vnmax, vntotal                          !セルの数、セルの数×2
            DO vn = 1,vnmax                                                             !セルの数だけループ
                write(n_unit,'(2(I12,2X))')  1, vn-1                                          !まずセルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_TYPES',vnmax
            DO vn = 1,vnmax                                                             !セルの数だけループ
                write(n_unit,'(I12)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行
            write(n_unit,'(A,I12)') 'CELL_DATA ',vnmax                                      !ここからセルのデータという合図、セルの数
            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !まずは飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1,vnmax                                                             !セルの数だけループ
                write(n_unit,'(F20.16)') radius(vn,1)*2.0d0                            !飛沫の直径
            END DO
            write(n_unit,'(A)') 'SCALARS Adhesion int'                                   !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1,vnmax                                                             !セルの数だけループ
                write(n_unit,'(I12)') adhesion(vn)                                            !次は飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO
            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1,vnmax                                                             !セルの数だけループ
                write(n_unit,'(3(F20.16,2X))') vel_drp(1,vn,1),vel_drp(2,vn,1),vel_drp(3,vn,1)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'WRITEOUT:', FN

    !=======飛沫データ（VTKファイル）の出力ここまで===========================

        !以下はCSVファイルの出力

        if(step==0) then !初期ステップならファイル新規作成
            open(newunit=n_unit, file=trim(path_out)//'/particle_'//trim(head_out)//'.csv', status='replace')
            print*,'REPLACE:particle_data.csv'

        else
            open(newunit=n_unit, file=trim(path_out)//'/particle_'//trim(head_out)//'.csv'&
                , action='write', status='old', position='append')

        end if
        
            write(n_unit,*) step*dt*L_chara/U_chara, ',', count(adhesion==0), ',', count(adhesion==2), ',', count(adhesion==3)
        close(n_unit)

        ! if(count(adhesion==0) <= 0) then !浮遊粒子がなくなれば計算終了
        !   print*,'all viruses terminated',N
        !   STOP
        ! end if

    end subroutine writeout

    !*************************************************************************************

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
        !-----calculation------------------------------------------
        TK = dble(T) + T0                                    ! 室温を絶対温度[K]に変換
        Es = Es0*exp((Lv/Rv)*(1.0d0/T0 - 1.0d0/TK))       ! 室温に置ける飽和蒸気圧

        norm = sqrt(sum(direction_g(:)**2))
        G(:) = G_dim * L_chara/(U_chara*U_chara) / norm * direction_g(:)
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
      
        if (radius(vn,1) <= rad_min(vn)) then  !半径が最小になったものを除く
                radius(vn,2) = rad_min(vn)
                return
        end if
    
        drdt1 = coeff / radius(vn,1)
        R1 = dt*drdt1
        R_approxi = radius(vn,1) + R1
        if(R_approxi <= 0.0d0) R_approxi = radius(vn,1)
        drdt2 = coeff / R_approxi
        R2 = dt*drdt2
        radius(vn,2) = radius(vn,1) + (R1+R2)*0.5d0
        
        if (radius(vn,2) < rad_min(vn)) radius(vn,2) = rad_min(vn)
      
        !*******************************************************************************************
    end subroutine evaporation
    !*******************************************************************************************

    !=====================================================================
    subroutine survival_check(step)
        integer,intent(in) :: step
        integer videal, vfloat, vn
            
        vfloat = count(adhesion == 0)
        videal = int(dble(vfloat)*(survival_rate(step) - survival_rate(step+1)))!videal:このステップで死滅すべき飛沫数
        vn = vnmax
    
        do while(videal > 0)
            if (adhesion(vn) == 0) then           !浮遊粒子からのみ除去する
                adhesion(vn) = -1
                crd_drp(:,vn,2) = -1.0d0
                vel_drp(:,vn,2) = 0.0d0

                vn = vn - 1
                videal = videal - 1
            else if (vn==0) then
                exit
            end if
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

    subroutine set_vn_trans(vnf)
        integer vn
        integer, intent(out) :: vnf
        vnf = 0
        do vn = 1, vnmax
            if(adhesion(vn) == 0) then
                vnf = vnf + 1
                vn_trans(vnf) = vn  !vn_trans(vnf):vnf(=1～vfloat)個目の浮遊飛沫粒子番号
            end if
        end do
    end subroutine set_vn_trans

    subroutine update_status(vn)
        integer, intent(in) :: vn
        !========= n step から n+1 step に入れ替え ===========================
        crd_drp(:,vn,1) = crd_drp(:,vn,2)
        radius(vn,1) = radius(vn,2)
        vel_drp(:,vn,1) = vel_drp(:,vn,2)
        !=====================================================================
    end subroutine update_status

    subroutine allocation_virus(num_virus)
        integer, intent(in) :: num_virus

        vnmax = num_virus

        allocate(crd_drp(3,num_virus,2),rad_min(num_virus),radius(num_virus,2), rad_ini(num_virus), vel_drp(3,num_virus,2))
        allocate(adhesion(num_virus),vn_trans(num_virus))

    end subroutine allocation_virus

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

end module virus_mod