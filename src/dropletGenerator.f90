module dropletGenerator_m
    use virusDroplet_m
    use dropletEquation_m
    implicit none
    private

    type placementBox
        double precision standardPoint(3), width(3)
    end type

    type, public :: DropletGenerator
        private

        type(DropletEquationSolver), pointer :: equation

        type(placementBox), allocatable :: pBox_array(:)

        double precision, allocatable :: radiusThreshold(:,:)

        integer :: generateRate = 0, generateMode = 0

        contains

        procedure, public :: generateDroplet
        procedure, public :: periodicGeneration => dropletPeriodicGeneration

        procedure set_dropletPlacementBox, set_dropletRadiusThreshold

        procedure calc_initialPosition
        procedure calc_initialRadius
        procedure set_virusDeadline

    end type

    integer, parameter :: nonActive = -99

    ! public BasicParameter, DropletEquationSolver, BasicParameter_, DropletEquationSolver_
    public DropletGroup, DropletGenerator_

    contains

    type(DropletGenerator) function DropletGenerator_( &
                    equation, radiusDistributionFile, positionDir, generationRate, generationMode)   !コンストラクタ
        type(DropletEquationSolver), target :: equation
        character(*), intent(in) :: radiusDistributionFile
        character(*), intent(in), optional :: positionDir
        integer, intent(in) :: generationRate, generationMode

        DropletGenerator_%equation => equation

        call DropletGenerator_%set_dropletRadiusThreshold(radiusDistributionFile)

        if(present(positionDir)) call DropletGenerator_%set_dropletPlacementBox(positionDir)

        DropletGenerator_%generateRate = generationRate
        DropletGenerator_%generateMode = generationMode
        select case(DropletGenerator_%generateMode)
        case(1)
            print*, "GenerationMode : [Generate]", DropletGenerator_%generateRate
        case(2)
            print*, "GenerationMode : [Activate]", DropletGenerator_%generateRate
        end select

    end function

    type(DropletGroup) function generateDroplet(self, num_droplet, nowTime, outputDir)
        class(DropletGenerator) self
        integer, intent(in) :: num_droplet
        double precision, intent(in) :: nowTime
        character(*), intent(in), optional :: outputDir

        if(num_droplet <= 0) return 

        allocate(generateDroplet%droplet(num_droplet))

        call self%calc_initialPosition(generateDroplet)

        if(present(outputDir)) then
            call self%calc_initialRadius(generateDroplet, outputDir)
        else
            call self%calc_initialRadius(generateDroplet)
        end if
        generateDroplet%droplet(:)%radius = generateDroplet%droplet(:)%initialRadius
        generateDroplet%droplet(:)%radius_min = self%equation%get_minimumRadius( &
                                                        generateDroplet%droplet(:)%initialRadius ) !最小半径の計算

        call self%set_virusDeadline(generateDroplet, nowTime)

        if(self%generateMode==2) generateDroplet%droplet(:)%status = nonActive

    end function

    subroutine calc_initialRadius(self, dGroup, outputDir)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        character(*), intent(in), optional :: outputDir
        integer vn, i, num_drop
        integer, allocatable :: rad_cnt(:)
        double precision random_rad, radius_dim(size(dGroup%droplet))

        num_drop = size(dGroup%droplet)

        allocate(rad_cnt(size(self%radiusThreshold, dim=2)), source=0)

        do vn = 1, num_drop                       !飛沫半径の分布を乱数によって与える

            call random_number(random_rad)

            do i = 1, size(self%radiusThreshold, dim=2)
                if(random_rad < self%radiusThreshold(2, i)) then
                    radius_dim(vn) = self%radiusThreshold(1, i)
                    rad_cnt(i) = rad_cnt(i) + 1
                    exit
                end if
            end do

        end do

        dGroup%droplet(:)%initialRadius = radius_dim(:) / self%equation%repValue('length') !初期飛沫半径のセットおよび無次元化

        if (sum(rad_cnt) /= num_drop) then
            print*, 'random_rad_ERROR', sum(rad_cnt), num_drop
            stop
        end if

        if(present(outputDir)) then
            block
                integer n_unit

                open(newunit=n_unit, file=outputDir//'/initialRadiusDistributon.csv', status='replace')
                    do i = 1, size(rad_cnt)
                        print '(4X, A, E10.3, A, I10)', 'rad_cnt(', self%radiusThreshold(1, i), ' ) =', rad_cnt(i)
                        write(n_unit, '(*(g0:,","))') self%radiusThreshold(1, i), rad_cnt(i)
                    end do
                close(n_unit)

            end block

        end if

    end subroutine

    subroutine calc_initialPosition(self, dGroup)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        integer kx,ky,kz, num_perEdge, num_perBox, k, k_end, cnt
        integer i_box, num_box, num_drop
        double precision :: standard(3), delta(3), width(3), randble(3)
        
        num_drop = size(dGroup%droplet)

        if(.not.allocated(self%pBox_array)) then
            print*, 'ERROR : InitialPositionBox is not Set.'
            stop
        end if
        num_box = size(self%pBox_array)

        num_perBox = num_drop / num_box
        
        print*, 'calc_initialPosition'

        num_perEdge = 1
        do while(num_box*((num_perEdge+1)**3) < num_drop)
            num_perEdge = num_perEdge + 1    !配置帯一辺当たりの飛沫数
        end do

        k = 1
        cnt = 1
        do i_box = 1, num_box
            width(:) = self%pBox_array(i_box)%width(:)
            standard(:) = self%pBox_array(i_box)%standardPoint(:)
  
            if(num_perEdge >= 2) then

                delta(:) = width(:) / dble(num_perEdge - 1)

                do kx = 1, num_perEdge

                    do ky = 1, num_perEdge

                        do kz = 1, num_perEdge

                            dGroup%droplet(k)%position(1) = standard(1) + delta(1)*dble(kx - 1)
                            dGroup%droplet(k)%position(2) = standard(2) + delta(2)*dble(ky - 1)
                            dGroup%droplet(k)%position(3) = standard(3) + delta(3)*dble(kz - 1)
                            k = k + 1
                            
                        end do

                    end do

                end do

            end if

            k_end = i_box * num_perBox
            if(i_box == num_box) k_end = num_drop

            do while(k <= k_end)
                call random_number(randble(:))
                dGroup%droplet(k)%position(:) = standard(:) + width(:)*randble(:)
                k = k + 1
            end do

            print*, 'BOX', i_box, 'has', k - cnt, 'droplets.'

            cnt = k

        end do

    end subroutine

    subroutine set_virusDeadline(self, dGroup, nowTime)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        double precision, intent(in) :: nowTime
        double precision randble(size(dGroup%droplet))

        call random_number(randble)

        dGroup%droplet(:)%deadline = self%equation%virusDeadline(randble(:)) + nowTime

    end subroutine

    subroutine set_dropletPlacementBox(self, positionDir)
        use simpleFile_reader
        use filename_mod
        class(DropletGenerator) self
        character(*), intent(in) :: positionDir
        integer i_box, num_box
        double precision, allocatable :: position_mat(:,:)
        character(:), allocatable :: fname
        logical existance
        
        fname = positionDir//'/'//IniPositionFName

        inquire(file=fname, exist=existance)

        if(.not.existance) then
            print*, '**Warning** '//fname//' is not found!'
            return
        end if

        call read_CSV(fname, position_mat)

        num_box = size(position_mat, dim=2)

        if(allocated(self%pBox_array)) deallocate(self%pBox_array)
        allocate(self%pBox_array(num_box))

        do i_box = 1, num_box
            self%pBox_array(i_box)%width(:) = position_mat(4:6, i_box)
            self%pBox_array(i_box)%standardPoint(:) = position_mat(1:3, i_box) - 0.5d0*self%pBox_array(i_box)%width(:)
        end do

    end subroutine

    subroutine set_dropletRadiusThreshold(self, radiusDistributionFilename)
        use simpleFile_reader
        class(DropletGenerator) self
        character(*), intent(in) :: radiusDistributionFilename

        if(allocated(self%radiusThreshold)) return

        call read_CSV('data/'//radiusDistributionFilename, self%radiusThreshold)

        self%radiusThreshold(1,:) = self%radiusThreshold(1,:) * 1.d-6   !マイクロメートル換算

    end subroutine

    function initialRadiusDistribution(self, dGroup) result(iniRadDis)
        class(DropletGenerator) self
        type(DropletGroup), intent(in) :: dGroup
        integer i, j, num_threshold
        real, allocatable :: iniRadDis(:,:)
        double precision threshold

        ! call set_dropletRadiusThreshold

        num_threshold = size(self%radiusThreshold, dim=2)
        allocate(iniRadDis(2,num_threshold))
        iniRadDis(1,:) = real(self%radiusThreshold(1,:))
        iniRadDis(2,:) = 0.0
        drop:do i = 1, size(dGroup%droplet)
            radius:do j = 1, num_threshold
                if(j < num_threshold) then
                    threshold = (self%radiusThreshold(1,j) + self%radiusThreshold(1,j+1))*0.5d0
                    if(dGroup%droplet(i)%initialRadius < threshold) then
                        iniRadDis(2,j) = iniRadDis(2,j) + 1.0
                        exit radius
                    end if
                else
                    iniRadDis(2,num_threshold) = iniRadDis(2,num_threshold) + 1.0

                end if
            end do radius
        end do drop

    end function

    subroutine dropletPeriodicGeneration(self, dGroup, nowTime)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        double precision, intent(in) :: nowTime
        integer num_generated, required_generation, num_nowGenerate

        required_generation = int(dble(self%generateRate)*nowTime)  !この時刻までに生成されているべき数

        select case(self%generateMode)
        case(1)
            num_generated = size(dGroup%droplet)  !今までに生成された数
            num_nowGenerate = required_generation - num_generated  !今このステップで生成されるべき数
            if(num_nowGenerate >= 1) then 
                dGroup = self%generateDroplet(num_nowGenerate, nowTime)
                dGroup%droplet = [dGroup%droplet, dGroup%droplet]
            end if

        case(2)
            num_generated = count(dGroup%droplet(:)%status /= nonActive)  !今までに生成された数
            if(num_generated == size(dGroup%droplet)) return        !生成され尽くした場合リターン

            num_nowGenerate = required_generation - num_generated !今このステップで生成されるべき数

            print*, num_generated, required_generation, num_nowGenerate
            
            block
                integer generateEnd, nonActive_perBox, i_box, num_box, generate_perBox
                integer, allocatable :: nonActiveID_array(:)

                num_box = size(self%pBox_array)

                if(num_nowGenerate >= num_box) then
                    nonActiveID_array = dGroup%IDinState(nonActive)

                    if(required_generation < size(dGroup%droplet)) then
                        nonActive_perBox = size(nonActiveID_array) / num_box
                        generate_perBox = num_nowGenerate / num_box

                        do i_box = 1, num_box
                            generateEnd = min(nonActive_perBox*(i_box-1) + generate_perBox, nonActive_perBox*i_box)

                            dGroup%droplet(nonActiveID_array(nonActive_perBox*(i_box-1) +1 : generateEnd))%status = 0

                        end do

                    else  !この時刻までに生成されているべき数が総飛沫数未満でない　＝＞　全て生成されるべき
                        dGroup%droplet(nonActiveID_array(:))%status = 0

                    end if
                
                end if
            end block

        end select

        ! print*, TimeOnSimu(), num_generated
    
    end subroutine


end module dropletGenerator_m