module equation_mod
    implicit none

    double precision, private :: delta_t !無次元時間間隔
    double precision, private :: L_represent, U_represent, Re

    double precision, parameter :: Rho_represent = 1.205d0          ! 空気の密度[kg/m3]
    double precision, parameter :: Mu_represent = 1.822d-5          ! 空気の粘性係数[kg/m3]
    double precision, parameter :: Rho_d = 0.99822d0*1.0d3          ! 飛沫（水）の密度[kg/m3]

    double precision, private :: coeff_drdt !半径変化率の無次元係数
    double precision, private :: gumma      !密度比（空気密度 / 飛沫(水)密度）
    double precision, private :: G(3)      !無次元重力加速度

    contains

    subroutine set_basical_variables(dt, L, U)
        double precision, intent(in) :: dt, L, U
        delta_t = dt
        L_represent = L
        U_represent = U

        gumma = Rho_represent / Rho_d     !密度比:    空気密度 / 飛沫(水)密度
        print*, 'gumma =', gumma

        Re = U_represent*L_represent*Rho_represent / Mu_represent

    end subroutine set_basical_variables

    subroutine set_gravity_acceleration(direction_g)
        double precision, intent(in) :: direction_g(3)
        double precision, parameter :: G_dim = 9.806650d0                          ! 重力加速度[m/s2]
        double precision norm

        norm = norm2(direction_g(:))
        G(:) = G_dim * L_represent/(U_represent*U_represent) / norm * direction_g(:)    !無次元重力加速度
        print*, 'Dimensionless Acceleration of Gravity :'
        print*, G(:)

    end subroutine set_gravity_acceleration

    !*************************************************************************************
    subroutine set_coeff_drdt(T, RH)
        !=====================================================================================
        real, intent(in) :: T, RH   !温度[℃]、相対湿度[%]
        double precision Es, TK
        double precision, parameter :: Rv = 461.51d0                           ! 水蒸気の気体定数[J/(kg.K)]
        double precision, parameter :: T0 = 273.15d0                               ! [K]
        double precision, parameter :: D = 0.2564d0*1.0d-4           ! 水蒸気の拡散定数[m2/s]
        double precision, parameter :: Lv = 2.451d0*1.0d6                  ! 水の蒸発潜熱[J/kg]
        double precision, parameter :: Es0 = 6.11d0*1.0d2                  ! 基準温度における飽和蒸気圧[Pa]
        !=====================================================================================  

        TK = dble(T) + T0                                    ! 室温を絶対温度[K]に変換
        Es = Es0*exp((Lv/Rv)*(1.0d0/T0 - 1.0d0/TK))       ! 室温に置ける飽和蒸気圧

        coeff_drdt = -D/(U_represent*L_represent) * (1.0d0 - dble(RH)/100.d0)*Es / (Rho_d*Rv*TK) ! dr/dt の無次元係数
        print*, 'coeff_drdt=', coeff_drdt

    end subroutine set_coeff_drdt
    !*******************************************************************************************

    function get_minimum_radius(initial_radius, RH) result(minimum_radius)
        use csv_reader
        implicit none
        double precision, intent(in) :: initial_radius(:)
        real, intent(in) :: RH
        double precision :: minimum_radius(size(initial_radius))
        double precision, allocatable :: rad_mat(:,:)
        integer i, i_max

        ! if(RH < 64) then
        !     minimum_radius(:) = initial_radius(:)*0.19d0
        ! else if(RH < 90) then
        !     minimum_radius(:) = initial_radius(:)*(0.073*exp(0.014*dble(RH)))
        ! else if(RH == 90) then
        !     minimum_radius(:) = initial_radius(:)*0.28d0
        ! else if(RH < 100) then
        !     minimum_radius(:) = initial_radius(:)*(0.0001*exp(0.0869*dble(RH)))
        ! else
        !     minimum_radius(:) = initial_radius(:)
        ! end if

        call read_CSV('data/minimum_radius.csv', rad_mat)
        i_max = size(rad_mat, dim=2)
        i = 1
        do while(rad_mat(1,i) < RH)
            i = i + 1
            if(i == i_max) exit
        end do
        print*, 'Dmin/D0 =', rad_mat(2,i), RH
        minimum_radius(:) = initial_radius(:) * rad_mat(2,i)

    end function get_minimum_radius

    !*******************************************************************************************
    double precision function evaporatin_eq(radius)
        double precision, intent(in) :: radius
        double precision drdt1,R1,R_approxi,drdt2,R2
      
        !========= 飛沫半径の変化の計算　(2次精度ルンゲクッタ（ホイン）) ===========================
    
        drdt1 = coeff_drdt / radius
        R1 = delta_t*drdt1

        R_approxi = radius + R1

        drdt2 = coeff_drdt / R_approxi
        R2 = delta_t*drdt2

        evaporatin_eq = radius + (R1+R2)*0.5d0
      
        !*******************************************************************************************
    end function evaporatin_eq
    !*******************************************************************************************

    !*******************************************************************************************
    function motion_eq(vel_d, vel_a, radius_d) result(vel_d_next)
        !*******************************************************************************************
        !=====================================================================================
        double precision, intent(in) :: vel_d(3), vel_a(3), radius_d
        double precision speed_r, Re_d, Cd, Coefficient, vel_d_next(3)
        !=====================================================================================

        speed_r = norm2(vel_a(:) - vel_d(:))    !相対速度の大きさ
        Re_d = (speed_r * 2.0d0*radius_d) * Re + 1.d-9  !ゼロ割回避のため、小さな値を足す

        Cd = (24.0d0/Re_d)*(1.0d0 + 0.15d0*(Re_d**0.687d0)) !抗力係数

        Coefficient = (3.0d0*Cd*gumma*speed_r)/(8.0d0*radius_d)

        vel_d_next(:) = ( vel_d(:) + ( G(:) + Coefficient*vel_a(:) )* delta_t ) &
                            / ( 1.0d0 + Coefficient * delta_t )

    end function motion_eq
    !----------------------------------------------------------------------------------

    !*******************************************************************************************
    function next_position(x1, v1, v2) result(x2)
        double precision, intent(in) :: x1(3), v1(3), v2(3)
        double precision x2(3)

        x2(:) = x1(:) + (v1(:) + v2(:))* 0.5d0 * delta_t

    end function next_position
    !----------------------------------------------------------------------------------
    
    double precision function survival_rate(step)
        integer, intent(in) :: step

        !このへんはインフルエンザのデータ（現在不使用）
        ! if(RH == 80)then  !　相対湿度80%の時使用
        !     survival_rate = 0.67d0*0.5102d0**(((L_represent/U_represent)*dt*dble(step-1))/3600.0d0)
        ! else if(RH == 50)then  !　相対湿度50%の時使用
        !     survival_rate = 0.84d0*0.5735d0**(((L_represent/U_represent)*dt*dble(step-1))/3600.0d0)
        ! else if(RH == 35)then  !　相対湿度35%の時使用
        !     survival_rate = 0.86d0*0.9240d0**(((L_represent/U_represent)*dt*dble(step-1))/3600.0d0)
        ! end if

        !新型コロナウイルス（1.1時間で半減）(論文によると、湿度30,60,90%のときのデータしかない)
        survival_rate = 0.999825d0**((L_represent/U_represent) * delta_t * dble(step-1))
    end function survival_rate

    double precision function representative_value(name)
        character(*), intent(in) :: name

        select case(name)
            case('Length')
                representative_value = L_represent
            case('Speed')
                representative_value = U_represent
            case('Time')
                representative_value = L_represent / U_represent

            case default
                representative_value = 0.0d0

        end select

    end function representative_value

    double precision function dimensional_time(step)
        integer, intent(in) :: step

        dimensional_time = step * delta_t * L_represent / U_represent

    end function dimensional_time

end module equation_mod