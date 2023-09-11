module dropletGenerator_m
    use virusDroplet_m
    use dropletEquation_m
    implicit none
    private

    type placementBox
        double precision center(3), width(3)
    end type

    type SequentialArray
        private
        integer index
        real, allocatable :: array(:)

        contains

        procedure set_SequentialArray
        procedure :: get_value => get_valueFromSequentialArray
        procedure :: get_valueArray => get_valueArrayFromSequentialArray

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
        procedure, public :: discharged_flag

        procedure set_dropletPlacementBox

        procedure calc_initialPosition

    end type

    ! public BasicParameter, DropletEquationSolver, BasicParameter_, DropletEquationSolver_
    public DropletGenerator_

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

    function generateDroplet(self, num_droplet, nowTime) result(droplets)
        class(DropletGenerator) self
        type(virusDroplet_t), allocatable :: droplets(:),transplant_droplets(:)
        integer, intent(in) :: num_droplet
        double precision, intent(in) :: nowTime
        double precision, allocatable :: initialRadius(:), deadline(:)

        integer :: k, kk, unit, num_droplet_transplant, cnt 

        num_droplet_transplant = 3248

        if(num_droplet <= 0) return 

        allocate(droplets(num_droplet))
        allocate(transplant_droplets(num_droplet_transplant))

        !call self%calc_initialPosition(droplets)
        !!!==連続計算での排出飛沫の初期位置を全部格納（0step to 3000000step）（2023/04/14〜）===!!!
        print*, '!!!!!!!!!!!!!!!!!!!!!calc_initialPosition!!!!!!!!!!!!!!!!!!!!!!!'

        open(newunit=unit,file='data/droplet_position.txt',status='old')
            do k = 1, num_droplet_transplant
                read(unit,'(3(f10.5))') transplant_droplets(k)%position(:)
            end do
        close(unit)

        cnt = 0
        do k = 1, num_droplet
            droplets(k)%position(:) = transplant_droplets(k - cnt*num_droplet_transplant)%position(:)
            if(mod(k,num_droplet_transplant) == 0) cnt = cnt + 1
        end do
        !!!=====================================================================================!!!

        initialRadius = self%initialRadiusArray%get_valueArray(num_droplet) * 1.d-6 !マイクロメートル換算
        initialRadius = initialRadius / self%equation%repValue('length') !無次元化
        !call set_initialRadius(droplets, initialRadius)
        !!!==連続計算での排出飛沫の飛沫半径を全部格納（0step to 3000000step）（2023/04/14〜）===!!!
        print*, '!!!!!!!!!!!!!!!!!!!!!set_initialRadius!!!!!!!!!!!!!!!!!!!!!!!'

        open(newunit=unit,file='data/droplet_diameter.txt',status='old')
            do k = 1, num_droplet_transplant
                read(unit,*) transplant_droplets(k)%initialRadius  !一度直径を読み込む
            end do
        close(unit)
        transplant_droplets(:)%initialRadius = transplant_droplets(:)%initialRadius *0.5d0      !直径を半径に変換
        transplant_droplets(:)%radius = transplant_droplets(:)%initialRadius

        cnt = 0
        do k = 1, num_droplet
            droplets(k)%initialRadius = transplant_droplets(k - cnt*num_droplet_transplant)%initialRadius
            droplets(k)%radius = transplant_droplets(k - cnt*num_droplet_transplant)%radius

            if(mod(k,num_droplet_transplant) == 0) cnt = cnt + 1
        end do

        !!!=====================================================================================!!!

        call set_radiusLowerLimit(droplets, self%equation%get_radiusLowerLimitRatio())

        deadline = self%deadlineArray%get_valueArray(num_droplet) / self%equation%repValue('time') !無次元化
        deadline = deadline + nowTime
        call set_virusDeadline(droplets, deadline)

        !!!==連続計算での排出飛沫の速度を全部格納（0step to 3000000step）（2023/04/14〜）===!!!
        print*, '!!!!!!!!!!!!!!!!!!!!!set_initialVelocity!!!!!!!!!!!!!!!!!!!!!!!'

        open(newunit=unit,file='data/droplet_velocity.txt',status='old')
            do k = 1, num_droplet_transplant
                read(unit,'(3(f10.5))') transplant_droplets(k)%velocity(:)
            end do
        close(unit)

        cnt = 0
        do k = 1, num_droplet
            droplets(k)%velocity(:) = transplant_droplets(k - cnt*num_droplet_transplant)%velocity(:)
            if(mod(k,num_droplet_transplant) == 0) cnt = cnt + 1
        end do
        !!!=====================================================================================!!!

        !!!==連続計算での排出飛沫の排出ステップを全部格納（0step to 3000000step）（2023/04/14〜）===!!!
        print*, '!!!!!!!!!!!!!!!!!!!!!set_discharged_step!!!!!!!!!!!!!!!!!!!!!!!'

        open(newunit=unit,file='data/discharged_step.txt',status='old')
            do k = 1, num_droplet_transplant
                read(unit,'(I7)') transplant_droplets(k)%discharged_step
            end do
        close(unit)

        cnt = 0
        do k = 1, num_droplet
            droplets(k)%discharged_step = transplant_droplets(k - cnt*num_droplet_transplant)%discharged_step + cnt*600 !移植の関係上100の倍数にしておく。1/10stepになるので注意!
            if(mod(k,num_droplet_transplant) == 0) cnt = cnt + 1
        end do
        !!!=====================================================================================!!!

        !if(self%generateRate > 0) 
        call set_dropletStatus(droplets, 'nonActive')       !移植するタイミングを指定するために一旦nonActiveにする

    end function

    subroutine calc_initialPosition(self, droplets)
        class(DropletGenerator), intent(in) :: self
        type(virusDroplet_t), intent(inout) ::  droplets(:)
        integer kx,ky,kz, num_perEdge, num_perBox, k, k_end, cnt
        integer i_box, num_box, num_drop
        double precision :: standard(3), delta(3), width(3) !randble(3)  

        if(.not.allocated(self%pBox_array)) then
            print*, 'ERROR : InitialPositionBox is not Set.'
            error stop
        end if
        num_box = size(self%pBox_array)

        num_drop = size(droplets)
        num_perBox = num_drop / num_box
        

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

                            droplets(k)%position(1) = standard(1) + delta(1)*dble(kx - 1)
                            droplets(k)%position(2) = standard(2) + delta(2)*dble(ky - 1)
                            droplets(k)%position(3) = standard(3) + delta(3)*dble(kz - 1)
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
                        ! droplets(k)%position(:) = standard(:) + width(:)*randble(:)
                        do d = 1, d_max
                            droplets(k)%position(:) = self%pBox_array(i_box)%center(:) + width(:)*dble(direction(:,d))
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

    subroutine set_dropletPlacementBox(self, positionDir)
        use simpleFile_reader
        use filename_m, only : IniPositionFName => InitialPositionFileName
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

    subroutine dropletPeriodicGeneration(self, droplets, nowTime, stat)
        class(DropletGenerator) self
        type(virusDroplet_t), intent(inout) :: droplets(:)
        double precision, intent(in) :: nowTime
        integer num_generated, required_generation, num_nowGenerate
        logical, intent(out) :: stat

        stat = .false.
        if(self%generateRate == 0) return
        
        required_generation = int(dble(self%generateRate)*nowTime * self%equation%repValue('time'))  !このステップ終了までに生成されているべき数

        num_generated = dropletCounter(droplets, 'total') - dropletCounter(droplets, 'nonActive')  !今までに生成された数
        if(num_generated == size(droplets)) return        !生成され尽くした場合リターン

        num_nowGenerate = required_generation - num_generated !今このステップで生成されるべき数

        ! print*, num_generated, required_generation, num_nowGenerate
        
        block
            integer generateEnd, nonActive_perBox, i_box, num_box, generate_perBox
            integer, allocatable :: nonActiveID_array(:)

            num_box = size(self%pBox_array)

            if(num_nowGenerate >= num_box) then
                nonActiveID_array = dropletIDinState(droplets, 'nonActive')

                if(required_generation < size(droplets)) then
                    nonActive_perBox = size(nonActiveID_array) / num_box
                    generate_perBox = num_nowGenerate / num_box

                    do i_box = 1, num_box
                        generateEnd = min(nonActive_perBox*(i_box-1) + generate_perBox, nonActive_perBox*i_box)
                        call set_dropletStatus(droplets, 'floating', nonActiveID_array(nonActive_perBox*(i_box-1) +1 : generateEnd))
                    end do

                else  !この時刻までに生成されているべき数が総飛沫数未満でない　＝＞　全て生成されるべき
                    call set_dropletStatus(droplets, 'floating', nonActiveID_array(:))

                end if

                stat = .true.
            
            end if
        end block

        ! print*, TimeOnSimu(), num_generated
    
    end subroutine

    subroutine discharged_flag(self, droplets, step, stat)
        class(DropletGenerator) :: self
        type(virusDroplet_t), intent(inout) :: droplets(:)
        integer, intent(in) :: step 
        integer, allocatable :: discharged_IDarray(:)
        integer :: vn, num_droplet,cnt
        logical, intent(out) :: stat
        logical, allocatable :: discharge_ID(:)

        stat = .false.
        if(allocated(discharged_IDarray)) deallocate(discharged_IDarray)

        num_droplet = size(droplets) 
        if(.not.allocated(discharge_ID)) allocate(discharge_ID(num_droplet))
        discharge_ID(:) = .false. 

        do vn = 1, num_droplet
            if(step == droplets(vn)%discharged_step) then 
                discharge_ID(vn) = .true.
            end if 
        end do   

        allocate(discharged_IDarray(count(discharge_ID)))
        cnt = 0 
        do vn = 1, num_droplet
            if(discharge_ID(vn)) then 
                cnt = cnt + 1 
                discharged_IDarray(cnt) = vn 
            end if 
        end do 

        call set_dropletStatus(droplets, 'floating', discharged_IDarray(:)) 

        stat = .true.

    end subroutine
        
    subroutine set_SequentialArray(self, filename)
        use array_m
        class(SequentialArray) self
        character(*), intent(in) :: filename

        call read_1dArray_real(filename, self%array)

        self%index = 1

    end subroutine

    real function get_valueFromSequentialArray(self)
        class(SequentialArray) self

        get_valueFromSequentialArray = self%array(self%index)

        self%index = self%index + 1
        if(self%index > size(self%array)) self%index = 1

    end function

    function get_valueArrayFromSequentialArray(self, arraySize) result(array)
        class(SequentialArray) self
        integer, intent(in) :: arraySize
        real, allocatable :: array(:)
        integer i

        allocate(array(arraySize))

        do i = 1, arraySize
            array(i) = self%get_value()
        end do

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