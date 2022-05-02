module dropletGenerator_m
    use virusDroplet_m
    use dropletEquation_m
    implicit none
    private

    type placementBox
        double precision center(3), width(3)
    end type

    type SequentialArray
        integer index
        real, allocatable :: array(:)

        contains

        procedure set_SequentialArray
        procedure :: get_value => get_valueFromSequentialArray

    end type

    type, public :: DropletGenerator
        private

        type(DropletEquationSolver), pointer :: equation

        type(placementBox), allocatable :: pBox_array(:)

        ! double precision, allocatable :: radiusThreshold(:,:)
        type(SequentialArray) initialRadiusArray, deadlineArray

        integer :: generateRate = 0

        contains

        procedure, public :: generateDroplet
        procedure, public :: periodicGeneration => dropletPeriodicGeneration

        procedure set_dropletPlacementBox

        procedure calc_initialPosition
        procedure set_initialRadius
        procedure set_virusDeadline

    end type

    integer, parameter :: nonActive = -99

    ! public BasicParameter, DropletEquationSolver, BasicParameter_, DropletEquationSolver_
    public DropletGroup, DropletGenerator_

    contains

    type(DropletGenerator) function DropletGenerator_( &
                    equation, radiusDistributionFile, positionDir, generationRate)   !コンストラクタ
        type(DropletEquationSolver), target :: equation
        character(*), intent(in) :: radiusDistributionFile
        character(*), intent(in) :: positionDir
        integer, intent(in) :: generationRate

        DropletGenerator_%equation => equation

        call DropletGenerator_%initialRadiusArray%set_SequentialArray('data/'//radiusDistributionFile)

        call DropletGenerator_%deadlineArray%set_SequentialArray('data/deadline.txt')

        call DropletGenerator_%set_dropletPlacementBox(positionDir)

        DropletGenerator_%generateRate = generationRate

    end function

    type(DropletGroup) function generateDroplet(self, num_droplet, nowTime)
        class(DropletGenerator) self
        integer, intent(in) :: num_droplet
        double precision, intent(in) :: nowTime

        if(num_droplet <= 0) return 

        allocate(generateDroplet%droplet(num_droplet))

        call self%calc_initialPosition(generateDroplet)

        call self%set_initialRadius(generateDroplet)

        generateDroplet%droplet(:)%radius = generateDroplet%droplet(:)%initialRadius
        generateDroplet%droplet(:)%radius_min = self%equation%get_minimumRadius( &
                                                        generateDroplet%droplet(:)%initialRadius ) !最小半径の計算

        call self%set_virusDeadline(generateDroplet, nowTime)

        if(self%generateRate > 0) generateDroplet%droplet(:)%status = nonActive

    end function

    subroutine set_initialRadius(self, dGroup)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        double precision initialRadius
        integer i

        do i = 1, size(dGroup%droplet)
            initialRadius = self%initialRadiusArray%get_value() * 1.d-6 !マイクロメートル換算
            dGroup%droplet(i)%initialRadius = initialRadius / self%equation%repValue('length') !無次元化
        end do

        ! block
        !     real, allocatable :: category(:)
        !     integer, allocatable :: frequency(:)
        !     call check_category(real(dGroup%droplet(:)%initialRadius), category, frequency)
        !     print*, category
        !     print*, frequency
        ! end block

    end subroutine

    subroutine calc_initialPosition(self, dGroup)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        integer kx,ky,kz, num_perEdge, num_perBox, k, k_end, cnt
        integer i_box, num_box, num_drop
        double precision :: standard(3), delta(3), width(3)!, randble(3)
        
        num_drop = size(dGroup%droplet)

        if(.not.allocated(self%pBox_array)) then
            print*, 'ERROR : InitialPositionBox is not Set.'
            stop
        end if
        num_box = size(self%pBox_array)

        num_perBox = num_drop / num_box
        
        ! print*, 'calc_initialPosition'

        num_perEdge = 1
        do while(num_box*((num_perEdge+1)**3) < num_drop)
            num_perEdge = num_perEdge + 1    !配置帯一辺当たりの飛沫数
        end do

        k = 1
        cnt = 1
        box:do i_box = 1, num_box
            k_end = i_box * num_perBox
            if(i_box == num_box) k_end = num_drop

            width(:) = self%pBox_array(i_box)%width(:)
            standard(:) = self%pBox_array(i_box)%center(:) - 0.5d0*self%pBox_array(i_box)%width(:)
  
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

            if(k <= k_end) then

                block
                    integer d, d_max
                    integer :: direction(3,6) = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], shape(direction))

                    width(:) = 0.d0
                    d_max = 1
                    placement:do
                        ! call random_number(randble(:))
                        ! dGroup%droplet(k)%position(:) = standard(:) + width(:)*randble(:)
                        do d = 1, d_max
                            dGroup%droplet(k)%position(:) = self%pBox_array(i_box)%center(:) + width(:)*dble(direction(:,d))
                            k = k + 1
                            if(k > k_end) exit placement
                        end do
                        width(:) = (width(:) + 0.5d0*self%pBox_array(i_box)%width(:)) * 0.5d0
                        d_max = 6
                    end do placement

                end block

            end if

            ! print*, 'BOX', i_box, 'has', k - cnt, 'droplets.'

            cnt = k

        end do box

    end subroutine

    subroutine set_virusDeadline(self, dGroup, nowTime)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        double precision, intent(in) :: nowTime
        double precision deadline
        integer i

        do i = 1, size(dGroup%droplet)
            deadline = self%deadlineArray%get_value() / self%equation%repValue('time') !無次元化
            dGroup%droplet(i)%deadline = deadline + nowTime
        end do

    end subroutine

    subroutine set_dropletPlacementBox(self, positionDir)
        use simpleFile_reader
        use filename_mod, only : IniPositionFName
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

        allocate(self%pBox_array(num_box))

        do i_box = 1, num_box
            self%pBox_array(i_box)%center(:) = position_mat(1:3, i_box)
            self%pBox_array(i_box)%width(:) = position_mat(4:6, i_box)
        end do

    end subroutine

    subroutine dropletPeriodicGeneration(self, dGroup, nowTime, stat)
        class(DropletGenerator) self
        type(DropletGroup) dGroup
        double precision, intent(in) :: nowTime
        integer num_generated, required_generation, num_nowGenerate
        logical, intent(out) :: stat

        stat = .false.
        required_generation = int(dble(self%generateRate)*nowTime * self%equation%repValue('time'))  !このステップ終了までに生成されているべき数

        num_generated = count(dGroup%droplet(:)%status /= nonActive)  !今までに生成された数
        if(num_generated == size(dGroup%droplet)) return        !生成され尽くした場合リターン

        num_nowGenerate = required_generation - num_generated !今このステップで生成されるべき数

        ! print*, num_generated, required_generation, num_nowGenerate
        
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

                stat = .true.
            
            end if
        end block

        ! print*, TimeOnSimu(), num_generated
    
    end subroutine

    subroutine set_SequentialArray(self, filename)
        use simpleFile_reader
        class(SequentialArray) self
        character(*), intent(in) :: filename

        call read_array_real(filename, self%array)

        self%index = 1

    end subroutine

    real function get_valueFromSequentialArray(self)
        class(SequentialArray) self

        get_valueFromSequentialArray = self%array(self%index)

        self%index = self%index + 1
        if(self%index > size(self%array)) self%index = 1

    end function

    subroutine check_category(array, categories, frequency)
        real, intent(in) :: array(:)
        real category_array(size(array))
        integer frequency_array(size(array))
        integer i, i_ctgry, num_category
        logical hit
        real, allocatable, intent(out) :: categories(:)
        integer, allocatable, intent(out) :: frequency(:)

        category_array = -1.e20

        category_array(1) = array(1)
        frequency_array(1) = 1
        num_category = 1

        do i = 2, size(array)
            hit = .false.
            do i_ctgry = 1, num_category
                if(array(i) == category_array(i_ctgry)) then
                    frequency_array(i_ctgry) = frequency_array(i_ctgry) + 1
                    hit = .true.
                    exit
                end if
            end do
            if(.not.hit) then
                num_category = num_category + 1
                category_array(num_category) = array(i)
                frequency_array(num_category) = 1
            end if
        end do

        categories = category_array(:num_category)
        frequency = frequency_array(:num_category)

    end subroutine

end module dropletGenerator_m