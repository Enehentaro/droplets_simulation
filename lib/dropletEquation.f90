module dropletEquation_m
    implicit none
    private

    type, public :: BasicParameter
        private
        double precision dt !無次元時間間隔
        double precision L, U, Re

        ! integer, public, target :: timeStep = 0   !構造体の要素はtarget属性にできないぽい

        contains
        procedure :: repValue => representativeValue
        procedure TimeStep2RealTime
    end type

    double precision, parameter :: Rho = 1.205d0          ! 空気の密度[kg/m3]
    double precision, parameter :: Mu = 1.822d-5          ! 空気の粘性係数[Pa・sec]
    double precision, parameter :: Rho_d = 0.99822d3          ! 飛沫（水）の密度[kg/m3]
    double precision, parameter :: gamma = Rho / Rho_d      !密度比（空気密度 / 飛沫(水)密度）

    type, public, extends(BasicParameter) :: DropletEquationSolver
        private
        double precision coeff_drdt !半径変化率の無次元係数
        double precision G(3)      !無次元重力加速度

        real T, RH
        double precision minimumRadiusRatio
        double precision, allocatable :: minimumRadiusMatrix(:,:)

        contains
        procedure set_gravity_acceleration, set_dropletEnvironment, dropletEnvironment
        procedure set_coeff_drdt, set_minimumRadiusRatio
        procedure next_position, next_velocity
        procedure get_radiusLowerLimitRatio

        procedure :: evaporationEq => evaporationEquation
        procedure solve_motionEquation

    end type

    public BasicParameter_, DropletEquationSolver_

    contains

    type(BasicParameter) function BasicParameter_(delta_t, L_represent, U_represent)
        double precision, intent(in) :: delta_t, L_represent, U_represent

        BasicParameter_%dt = delta_t
        BasicParameter_%L = L_represent
        BasicParameter_%U = U_represent

        print*, 'Delta_Time =', BasicParameter_%dt

        BasicParameter_%Re = BasicParameter_%U*BasicParameter_%L*Rho / Mu

    end function

    type(DropletEquationSolver) function DropletEquationSolver_(delta_t, L_represent, U_represent, &
                                                                direction_g, Temperature, RelativeHumidity)
        double precision, intent(in) :: delta_t, L_represent, U_represent
        double precision, intent(in) :: direction_g(3)
        real, intent(in) :: Temperature, RelativeHumidity

        DropletEquationSolver_%BasicParameter = BasicParameter_(delta_t, L_represent, U_represent)

        call DropletEquationSolver_%set_gravity_acceleration(direction_g)
        call DropletEquationSolver_%set_dropletEnvironment(Temperature, RelativeHumidity)

    end function

    subroutine set_gravity_acceleration(self, direction_g)
        use vector_m
        class(DropletEquationSolver) self
        double precision, intent(in) :: direction_g(3)
        double precision, parameter :: G_dim = 9.806650d0                          ! 重力加速度[m/s2]
        double precision norm

        norm = G_dim * self%L/(self%U*self%U)
        self%G(:) = norm * normalize_vector(direction_g(:))    !無次元重力加速度
        print*, 'Dimensionless Acceleration of Gravity :'
        print*, self%G(:)

    end subroutine

    subroutine set_dropletEnvironment(self, Temperature, RelativeHumidity)
        class(DropletEquationSolver) self
        real, intent(in) :: Temperature, RelativeHumidity

        self%T = Temperature
        self%RH = RelativeHumidity

        call self%set_coeff_drdt()          !温湿度依存の係数の設定
        call self%set_minimumRadiusRatio()

    end subroutine

    real function dropletEnvironment(self, name)
        class(DropletEquationSolver) self
        character(*), intent(in) :: name

        select case(name)
            case('Temperature')
                dropletEnvironment = self%T
            case('RelativeHumidity')
                dropletEnvironment = self%RH
            case default
                dropletEnvironment = -1.e20
        end select

    end function

    subroutine set_coeff_drdt(self)
        class(DropletEquationSolver) self
        !=====================================================================================
        double precision Es, TK
        double precision, parameter :: Rv = 461.51d0                           ! 水蒸気の気体定数[J/(kg.K)]
        double precision, parameter :: T0 = 273.15d0                               ! [K]
        double precision, parameter :: D = 0.2564d-4           ! 水蒸気の拡散定数[m2/s]
        double precision, parameter :: Lv = 2.451d6                  ! 水の蒸発潜熱[J/kg]
        double precision, parameter :: Es0 = 6.11d2                  ! 基準温度における飽和蒸気圧[Pa]
        !=====================================================================================  
        TK = dble(self%T) + T0                                    ! 室温を絶対温度[K]に変換
        Es = Es0*exp((Lv/Rv)*(1.0d0/T0 - 1.0d0/TK))       ! 室温における飽和蒸気圧

        self%coeff_drdt = -D/(self%U*self%L) * (1.0d0 - dble(self%RH)/100.d0)*Es / (Rho_d*Rv*TK) ! dr/dt の無次元係数
        print*, 'coeff_drdt=', self%coeff_drdt

    end subroutine

    subroutine set_minimumRadiusRatio(self)
        use simpleFile_reader
        class(DropletEquationSolver) self
        integer i, i_max

        call read_CSV('data/minimum_radius.csv', self%minimumRadiusMatrix)
        
        i_max = size(self%minimumRadiusMatrix, dim=2)
        i = 1
        do while(self%minimumRadiusMatrix(1,i) < self%RH)
            i = i + 1
            if(i == i_max) exit
        end do
        self%minimumRadiusRatio = self%minimumRadiusMatrix(2,i)

        print*, 'Dmin/D0 =', self%minimumRadiusRatio, self%RH

    end subroutine

    double precision function get_radiusLowerLimitRatio(self)
        class(DropletEquationSolver), intent(in) :: self

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

        get_radiusLowerLimitRatio = self%minimumRadiusRatio

    end function

    !蒸発方程式。半径変化量を返す。
    function evaporationEquation(self, radius) result(dr)
        class(DropletEquationSolver) self
        double precision, intent(in) :: radius
        double precision drdt1,dr1, drdt2,dr2, r_approxi, dr
        !========= 飛沫半径の変化の計算　(2次精度ルンゲクッタ（ホイン）) ===========================
    
        drdt1 = self%coeff_drdt / radius
        dr1 = drdt1 * self%dt

        r_approxi = radius + dr1

        drdt2 = self%coeff_drdt / r_approxi
        dr2 = drdt2 * self%dt

        dr = (dr1 + dr2)*0.5d0

    end function

    subroutine solve_motionEquation(self, X, V, Va, R)
        class(DropletEquationSolver) self
        double precision, intent(inout) :: X(3), V(3)
        double precision, intent(in) :: Va(3), R
        double precision V_now(3)

        V_now(:) = V(:)

        V(:) = self%next_velocity(V_now(:), Va(:), R)

        X(:) = self%next_position(X(:), V_now(:), V(:))
    
    end subroutine

    function next_velocity(self, vel_d, vel_a, radius_d) result(vel_d_next)
        class(DropletEquationSolver) self
        double precision, intent(in) :: vel_d(3), vel_a(3), radius_d
        double precision speed_r, Re_d, CD, C, vel_d_next(3)

        if(radius_d <= 0.d0) then
            print*, '**ZeroRadius ERROR**', radius_d
            error stop
        end if

        speed_r = norm2(vel_a(:) - vel_d(:))    !相対速度の大きさ
        Re_d = (speed_r * 2.0d0*radius_d) * self%Re

        CD = DragCoefficient(Re_d) !抗力係数

        C = (3.0d0*CD*gamma*speed_r)/(8.0d0*radius_d)

        vel_d_next(:) = ( vel_d(:) + ( self%G(:) + C*vel_a(:) )* self%dt ) / ( 1.0d0 + C*self%dt )

    end function

    function next_position(self, x1, v1, v2) result(x2)
        class(DropletEquationSolver) self
        double precision, intent(in) :: x1(3), v1(3), v2(3)
        double precision x2(3)

        x2(:) = x1(:) + (v1(:) + v2(:))* 0.5d0 * self%dt

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

    ! elemental double precision function virusDeadline(self, deathParameter)
    !     class(DropletEquationSolver), intent(in) :: self
    !     double precision, intent(in) :: deathParameter
    !     double precision, parameter :: halfLife = 3960.d0   !半減期 1.1 h ( = 3960 sec)
    !     double precision, parameter :: alpha = log(2.d0) / halfLife

    !     virusDeadline = - log(deathParameter) / alpha
    !     virusDeadline = virusDeadline * self%U/self%L !無次元化

    ! end function

    double precision function representativeValue(self, name)
        class(BasicParameter) self
        character(*), intent(in) :: name

        ! if(L*U <= 0.d0) print*, '**WARNING** ZeroRepresentativeValue :', L, U   !代表値がゼロなら警告

        select case(name)
            case('length')
                representativeValue = self%L
            case('speed')
                representativeValue = self%U
            case('time')
                representativeValue = self%L / self%U

            case default
                print*, 'RepresentativeValueERROR : ', name
                error stop

        end select

    end function

    double precision function TimeStep2RealTime(self, step, dimension)
        class(BasicParameter) self
        integer, intent(in) :: step
        logical, intent(in) :: dimension

        TimeStep2RealTime = step * self%dt
        
        if(dimension) TimeStep2RealTime = TimeStep2RealTime * self%repValue('time')

    end function

    ! double precision function deltaTime()

    !     deltaTime = dt
        
    ! end function

end module dropletEquation_m