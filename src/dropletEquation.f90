module dropletEquation_m
    implicit none
    private

    double precision dt !無次元時間間隔
    double precision L, U, Re

    double precision, parameter :: Rho = 1.205d0          ! 空気の密度[kg/m3]
    double precision, parameter :: Mu = 1.822d-5          ! 空気の粘性係数[Pa・sec]
    double precision, parameter :: Rho_d = 0.99822d3          ! 飛沫（水）の密度[kg/m3]
    double precision, parameter :: gumma = Rho / Rho_d      !密度比（空気密度 / 飛沫(水)密度）

    double precision coeff_drdt !半径変化率の無次元係数
    double precision G(3)      !無次元重力加速度

    real T, RH
    double precision minimumRadiusRatio
    double precision, allocatable :: minimumRadiusMatrix(:,:)

    public set_basicVariables_dropletEquation, set_gravity_acceleration, set_dropletEnvironment, dropletEnvironment
    public evaporatin_eq, solve_motionEquation, representativeValue, deltaTime
    public get_minimumRadius, virusDeadline

    contains

    subroutine set_basicVariables_dropletEquation(delta_t, L_represent, U_represent)
        double precision, intent(in) :: delta_t, L_represent, U_represent

        dt = delta_t
        L = L_represent
        U = U_represent

        print*, 'Delta_Time =', dt

        Re = U*L*Rho / Mu

    end subroutine

    subroutine set_gravity_acceleration(direction_g)
        use vector_m
        double precision, intent(in) :: direction_g(3)
        double precision, parameter :: G_dim = 9.806650d0                          ! 重力加速度[m/s2]
        double precision norm

        norm = G_dim * L/(U*U)
        G(:) = norm * normalize_vector(direction_g(:))    !無次元重力加速度
        print*, 'Dimensionless Acceleration of Gravity :'
        print*, G(:)

    end subroutine

    subroutine set_dropletEnvironment(Temperature, RelativeHumidity)
        real, intent(in) :: Temperature, RelativeHumidity

        T = Temperature
        RH = RelativeHumidity

        call set_coeff_drdt          !温湿度依存の係数の設定
        call set_minimumRadiusRatio

    end subroutine

    real function dropletEnvironment(name)
        character(*), intent(in) :: name

        select case(name)
            case('Temperature')
                dropletEnvironment = T
            case('RelativeHumidity')
                dropletEnvironment = RH
            case default
                dropletEnvironment = -1.e20
        end select

    end function

    subroutine set_coeff_drdt
        !=====================================================================================
        double precision Es, TK
        double precision, parameter :: Rv = 461.51d0                           ! 水蒸気の気体定数[J/(kg.K)]
        double precision, parameter :: T0 = 273.15d0                               ! [K]
        double precision, parameter :: D = 0.2564d-4           ! 水蒸気の拡散定数[m2/s]
        double precision, parameter :: Lv = 2.451d6                  ! 水の蒸発潜熱[J/kg]
        double precision, parameter :: Es0 = 6.11d2                  ! 基準温度における飽和蒸気圧[Pa]
        !=====================================================================================  
        TK = dble(T) + T0                                    ! 室温を絶対温度[K]に変換
        Es = Es0*exp((Lv/Rv)*(1.0d0/T0 - 1.0d0/TK))       ! 室温における飽和蒸気圧

        coeff_drdt = -D/(U*L) * (1.0d0 - dble(RH)/100.d0)*Es / (Rho_d*Rv*TK) ! dr/dt の無次元係数
        print*, 'coeff_drdt=', coeff_drdt

    end subroutine

    subroutine set_minimumRadiusRatio
        use simpleFile_reader
        integer i, i_max

        if(allocated(minimumRadiusMatrix)) return

        call read_CSV('data/minimum_radius.csv', minimumRadiusMatrix)
        
        i_max = size(minimumRadiusMatrix, dim=2)
        i = 1
        do while(minimumRadiusMatrix(1,i) < RH)
            i = i + 1
            if(i == i_max) exit
        end do
        minimumRadiusRatio = minimumRadiusMatrix(2,i)

        print*, 'Dmin/D0 =', minimumRadiusRatio, RH

    end subroutine

    elemental double precision function get_minimumRadius(initial_radius)
        double precision, intent(in) :: initial_radius

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

        get_minimumRadius = initial_radius * minimumRadiusRatio

    end function

    double precision function evaporatin_eq(radius)
        double precision, intent(in) :: radius
        double precision drdt1,dr1, drdt2,dr2, r_approxi
        !========= 飛沫半径の変化の計算　(2次精度ルンゲクッタ（ホイン）) ===========================
    
        drdt1 = coeff_drdt / radius
        dr1 = drdt1 * dt

        r_approxi = radius + dr1

        drdt2 = coeff_drdt / r_approxi
        dr2 = drdt2 * dt

        evaporatin_eq = radius + (dr1 + dr2)*0.5d0

    end function

    subroutine solve_motionEquation(X, V, Va, R)
        double precision, intent(inout) :: X(3), V(3)
        double precision, intent(in) :: Va(3), R
        double precision V_now(3)

        V_now(:) = V(:)

        V(:) = next_velocity(V_now(:), Va(:), R)

        X(:) = next_position(X(:), V_now(:), V(:))
    
    end subroutine

    function next_velocity(vel_d, vel_a, radius_d) result(vel_d_next)
        double precision, intent(in) :: vel_d(3), vel_a(3), radius_d
        double precision speed_r, Re_d, CD, C, vel_d_next(3)

        if(radius_d <= 0.d0) then
            print*, '**zeroRadius ERROR**', radius_d
            stop
        end if

        speed_r = norm2(vel_a(:) - vel_d(:))    !相対速度の大きさ
        Re_d = (speed_r * 2.0d0*radius_d) * Re

        CD = DragCoefficient(Re_d) !抗力係数

        C = (3.0d0*CD*gumma*speed_r)/(8.0d0*radius_d)

        vel_d_next(:) = ( vel_d(:) + ( G(:) + C*vel_a(:) )* dt ) / ( 1.0d0 + C*dt )

    end function

    function next_position(x1, v1, v2) result(x2)
        double precision, intent(in) :: x1(3), v1(3), v2(3)
        double precision x2(3)

        x2(:) = x1(:) + (v1(:) + v2(:))* 0.5d0 * dt

    end function

    double precision function DragCoefficient(Re_d)
        double precision, intent(in) :: Re_d
        double precision Re_d_

        Re_d_ = Re_d + 1.d-9  !ゼロ割回避のため、小さな値を足す

        DragCoefficient = (24.0d0/Re_d_)*(1.0d0 + 0.15d0*(Re_d_**0.687d0))

    end function
    
    ! double precision function survival_rate(step)
    !     integer, intent(in) :: step
    !     double precision time

    !     !このへんはインフルエンザのデータ（現在不使用）
    !     ! if(RH == 80)then  !　相対湿度80%の時使用
    !     !     survival_rate = 0.67d0*0.5102d0**(((L/U)*dt*dble(step-1))/3600.0d0)
    !     ! else if(RH == 50)then  !　相対湿度50%の時使用
    !     !     survival_rate = 0.84d0*0.5735d0**(((L/U)*dt*dble(step-1))/3600.0d0)
    !     ! else if(RH == 35)then  !　相対湿度35%の時使用
    !     !     survival_rate = 0.86d0*0.9240d0**(((L/U)*dt*dble(step-1))/3600.0d0)
    !     ! end if

    !     time = TimeOnSimu(step, dimension=.true.)
    !     !新型コロナウイルス（1.1時間で半減）
    !     survival_rate = 0.999825d0**(time)
    ! end function

    elemental double precision function virusDeadline(deathParameter)
        double precision, intent(in) :: deathParameter
        double precision, parameter :: halfLife = 3960.d0   !半減期 1.1 h ( = 3960 sec)
        double precision, parameter :: alpha = log(2.d0) / halfLife

        virusDeadline = - log(deathParameter) / alpha
        virusDeadline = virusDeadline * U/L !無次元化

    end function

    double precision function representativeValue(name)
        character(*), intent(in) :: name

        select case(name)
            case('length')
                representativeValue = L
            case('speed')
                representativeValue = U
            case('time')
                representativeValue = L / U

            case default
                representativeValue = -1.d20

        end select

    end function

    double precision function deltaTime()

        deltaTime = dt
        
    end function

end module dropletEquation_m