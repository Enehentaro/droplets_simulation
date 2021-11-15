module virusDroplet_m
    implicit none

    type virusDroplet_t
        double precision :: position(3), velocity(3)=0.d0
        double precision radius, radius_min, deathParam
        integer :: status=0
    end type virusDroplet_t

    type(virusDroplet_t), allocatable, private :: droplets_ini(:)

    integer, allocatable :: leaderID(:)

    integer, allocatable :: statusCSV(:)

    contains

    subroutine allocation_initialDroplets(num_drop)
        integer, intent(in) :: num_drop

        if(allocated(droplets_ini)) deallocate(droplets_ini)
        allocate(droplets_ini(num_drop))

    end subroutine allocation_initialDroplets

    subroutine calc_initial_radius
        use equation_mod
        use csv_reader
        integer vn, i, num_drop
        double precision, allocatable :: threshold(:,:)
        integer, allocatable :: rad_cnt(:)
        double precision :: radius_dim(size(droplets_ini))
        real(8) random_rad

        num_drop = size(droplets_ini)
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

        droplets_ini(:)%radius = radius_dim(:) / representative_value('Length') !初期飛沫半径のセットおよび無次元化

        if (sum(rad_cnt) /= num_drop) then
            print*, 'random_rad_ERROR', sum(rad_cnt), num_drop
            stop
        end if

        do i = 1, size(rad_cnt)
            print*,'rad_cnt(', threshold(1, i), ') =', rad_cnt(i)
        end do

    end subroutine calc_initial_radius

    subroutine calc_initial_position(dir)
        use csv_reader
        use filename_mod
        character(*), intent(in) :: dir
        integer kx,ky,kz, num_per_edge, num_per_box, m, k, k_end, cnt
        integer i_box, num_box, num_drop
        double precision :: standard(3), delta(3), width(3), randble(3)
        double precision, allocatable :: position_mat(:,:)
        
        num_drop = size(droplets_ini)

        call read_CSV(dir//'/'//IniPositionFName, position_mat)

        num_box = size(position_mat, dim=2)

        num_per_box = num_drop / num_box
        
        print*, 'calc_initial_position'

        m = 1
        do while(num_box*(m**3) < num_drop)
            m = m + 1
        end do
        num_per_edge = m - 1    !配置帯一辺当たりの飛沫数

        if(allocated(leaderID)) deallocate(leaderID)
        allocate(leaderID(num_box + 1))

        k = 1
        cnt = 1
        do i_box = 1, num_box
            leaderID(i_box) = k
  
            width(:) = position_mat(4:6, i_box)
            standard(:) = position_mat(1:3, i_box) - 0.5d0*width(:)
            delta(:) = width(:) / dble(num_per_edge - 1)

            do kx = 1, num_per_edge

                do ky = 1, num_per_edge

                    do kz = 1, num_per_edge

                        droplets_ini(k)%position(1) = standard(1) + delta(1)*dble(kx - 1)
                        droplets_ini(k)%position(2) = standard(2) + delta(2)*dble(ky - 1)
                        droplets_ini(k)%position(3) = standard(3) + delta(3)*dble(kz - 1)
                        k = k + 1
                        
                    end do

                end do

            end do

            k_end = i_box * num_per_box
            if(i_box == num_box) k_end = num_drop

            do while(k <= k_end)
                call random_number(randble(:))
                droplets_ini(k)%position(:) = standard(:) + width(:)*randble(:)
                k = k + 1
            end do

            print*, 'BOX', i_box, 'has', k - cnt, 'droplets. LeaderID:', leaderID(i_box)

            cnt = k

        end do

        leaderID(num_box + 1) = num_drop + 1

    end subroutine calc_initial_position

    subroutine set_deathParam
        integer i, num_drop
        double precision randble

        num_drop = size(droplets_ini)

        do i = 1, num_drop
            call random_number(randble)
            droplets_ini(i)%deathParam = randble
        end do

    end subroutine set_deathParam

    subroutine read_initialDistribution(dir)
        character(*), intent(in) :: dir

        droplets_ini = read_droplet_VTK(dir//'/InitialDistribution.vtk') !自動割付
        leaderID = [1, size(droplets_ini)+1]

    end subroutine read_initialDistribution

    subroutine calc_minimumRadius(RH)
        use equation_mod
        real, intent(in) :: RH

        droplets_ini(:)%radius_min = get_minimum_radius(droplets_ini(:)%radius, RH) !最小半径の計算

    end subroutine calc_minimumRadius

    function get_initialState_of_droplets() result(initialDroplets)
        type(virusDroplet_t), allocatable :: initialDroplets(:)

        initialDroplets = droplets_ini
            
    end function get_initialState_of_droplets

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

    function read_droplet_VTK(fname) result(droplets_read)
        implicit none
        character(*), intent(in) :: fname
        type(virusDroplet_t), allocatable :: droplets_read(:)
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
                read(n_unit, *) droplets_read(vn)%position(:)
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
                read(n_unit,'(I12)') droplets_read(vn)%status
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) diameter(vn)
            END DO

            read(n_unit,'()')
            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets_read(vn)%deathParam
            END DO

            read(n_unit,'()')
            DO vn = 1, num_drop
                read(n_unit, *) droplets_read(vn)%velocity(:)
            END DO
        close(n_unit)
    
        droplets_read(:)%radius = diameter(:) * 0.5d0
      
    end function read_droplet_VTK

    subroutine output_droplet_VTK(fname, droplets, initial)
        implicit none
        character(*), intent(in) :: fname
        type(virusDroplet_t), intent(in) :: droplets(:)
        logical, optional :: initial
        integer vn, n_unit, num_drop

        num_drop = size(droplets)
        !=======ここから飛沫データ（VTKファイル）の出力===========================
        open(newunit=n_unit, file=fname, status='replace')                                             !ここで出力ファイルを指定
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'                                !ファイルの始め4行は文字列（決まり文句）
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ', num_drop,' float'                              !節点の数
            DO vn = 1, num_drop                                                            !節点の数だけループ
                write(n_unit,'(3(f20.15,2X))') droplets(vn)%position(:)   !節点の座標（左から順にx,y,z）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12,2X,I12)') 'CELLS ', num_drop, num_drop*2                          !セルの数、セルの数×2
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(2(I12,2X))')  1, vn-1                                          !まずセルを構成する点の数（セル形状が点なので1）、その点のID
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12)') 'CELL_TYPES', num_drop
            DO vn = 1, num_drop                                                           !セルの数だけループ
                write(n_unit,'(I12)') 1                                                        !セルの形状（1は点であることを意味する）
            END DO
            write(n_unit,'()')                                                              !改行

            write(n_unit,'(A,I12)') 'CELL_DATA ', num_drop                                      !ここからセルのデータという合図、セルの数
            
            write(n_unit,'(A)') 'SCALARS Status int'                                   !飛沫の状態(0:浮遊、1:付着、2:回収)
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(I12)') droplets(vn)%status                                            !飛沫の状態(0:浮遊、1:付着、2:回収)
            END DO

            write(n_unit,'(A)') 'SCALARS Diameter float'                                  !飛沫の直径
            write(n_unit,'(A)') 'LOOKUP_TABLE default'
            DO vn = 1, num_drop                                                            !セルの数だけループ
                write(n_unit,'(f20.15)') droplets(vn)%radius*2.0d0                            !飛沫の直径
            END DO

            if(present(initial)) then
                if(initial) then
                    write(n_unit,'(A)') 'SCALARS DeathParam float'
                    write(n_unit,'(A)') 'LOOKUP_TABLE default'
                    DO vn = 1, num_drop                                                            !セルの数だけループ
                        write(n_unit,'(f20.15)') droplets(vn)%deathParam
                    END DO
                end if
            end if

            write(n_unit,'(A)') 'VECTORS Velocity float'                             !最後に飛沫の速度
            DO vn = 1, num_drop                                                             !セルの数だけループ
                write(n_unit,'(3(f20.15,2X))') droplets(vn)%velocity(:)               !飛沫の速度
            END DO
        close(n_unit)

        print*, 'writeOUT:', trim(fname)

        !=======飛沫データ（VTKファイル）の出力ここまで===========================
    end subroutine output_droplet_VTK

    subroutine output_droplet_CSV(fname, droplets, time, initial)
        character(*), intent(in) :: fname
        type(virusDroplet_t), intent(in) :: droplets(:)
        double precision, intent(in) :: time
        logical, intent(in) :: initial
        integer n_unit, L

        if(initial) then !初回ならファイル新規作成
            open(newunit=n_unit, file=fname, status='replace')

        else
            open(newunit=n_unit, file=fname, action='write', status='old', position='append')

        end if

        if(.not.allocated(statusCSV)) statusCSV = [0, 1, -1,-2]
        
        write(n_unit,'(*(g0:,","))') real(time), (count(droplets(:)%status==statusCSV(L)), L = 1, size(statusCSV))

        close(n_unit)

    end subroutine output_droplet_CSV

end module virusDroplet_m