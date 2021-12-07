!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use dropletMotionSimulation
    implicit none
    ! type(dropletGroup) dGroup
    ! double precision, save :: next_time = 0.d0
    type(dropletGroup) dGroup
    double precision, save :: next_time = 0.d0
    integer n_unit
    integer cnt1,cnt2,cnt3,cnt4,cnt5,cnt6,sum
    integer f,g,h,i,j,k
    character(50) fname
    double precision,allocatable :: left(:),right(:),w(:),x(:),y(:),z(:)

    ! if(mod(timeStep, 100) == 0) then    !100ステップおきに飛沫を10個発生させるサンプルコード
    !     dGroup = generate_dropletGroup(10)
    !     call mainDroplet%append(dGroup)
    ! end if

    ! if(TimeOnSimu(dimension=.true.) >= next_time) then    !0.5秒おきに飛沫を150個発生させるサンプルコード
    !     dGroup = generate_dropletGroup(150)
    !     call mainDroplet%append(dGroup)
    !     next_time = next_time + 0.5d0
    ! end if

    ! dGroup = mainDroplet%inBox([-1.d0,-1.d0,-1.d0], [1.d0,1.d0,1.d0])  !ボックス内の飛沫を取得
    ! print*, size(dGroup%droplet)
    
    if(TimeOnSimu(dimension=.true.) >= next_time) then    !0.5秒おきに飛沫を150個発生させるサンプルコード
        dGroup = generate_dropletGroup(10)
        call mainDroplet%append(dGroup)
        next_time = next_time + 0.1d0
    end if

    ! lung inlet =======================================
    dGroup = mainDroplet%inBox([-0.036d0,-0.4629d0,0.952d0], [-0.0329d0,-0.4551d0,0.9602d0])
    cnt1 = size(dGroup%droplet)

    allocate(left(cnt1))
    
    do f = 1, cnt1
        left(f) = dGroup%droplet(f)%initialRadius
    end do

    dGroup = mainDroplet%inBox([0.0075d0,-0.4631d0,0.957d0], [0.0155d0,-0.4549d0,0.9608d0])
    cnt2 = size(dGroup%droplet)

    allocate(right(cnt2))
    
    do g = 1, cnt2
        right(g) = dGroup%droplet(g)%initialRadius
    end do

    ! =====================================================

    dGroup = mainDroplet%inBox([-0.036d0,-0.468d0,0.952d0], [0.016d0,-0.447d0,0.982d0])  !ボックス内の飛沫を取得
    cnt3 = size(dGroup%droplet)

    allocate(w(cnt3))
    
    do h = 1, cnt3
        w(h) = dGroup%droplet(h)%initialRadius
    end do

    dGroup = mainDroplet%inBox([-0.020d0,-0.468d0,0.982d0], [0.031d0,-0.433d0,1.104d0]) 
    cnt4 = size(dGroup%droplet)

    allocate(x(cnt4))
    
    do i = 1, cnt4
        x(i) = dGroup%droplet(i)%initialRadius
    end do
    
    dGroup = mainDroplet%inBox([-0.014d0,-0.468d0,1.104d0], [0.014d0,-0.3754d0,1.1289d0]) 
    cnt5 = size(dGroup%droplet)

    allocate(y(cnt5))
    
    do j = 1, cnt5
        y(j) = dGroup%droplet(j)%initialRadius
    end do

    dGroup = mainDroplet%inBox([-0.011d0,-0.468d0,1.1289d0], [0.011d0,-0.3601d0,1.1632d0]) 
    cnt6 = size(dGroup%droplet)

    allocate(z(cnt6))
    
    do k = 1, cnt6
        z(k) = dGroup%droplet(k)%initialRadius
    end do

    if (mod(timeStep, outputInterval) == 0) then
        sum = cnt3 + cnt4 + cnt5 + cnt6
        print*, 'cnt1 = ', cnt1
        print*, 'cnt2 = ', cnt2
        print*, 'cnt3 = ', cnt3
        print*, 'cnt4 = ', cnt4
        print*, 'cnt5 = ', cnt5
        print*, 'cnt6 = ', cnt6
        print*, 'sum  = ', sum
    end if

    if (timeStep == n_end) then

    fname = 'output14.txt'
    open(newunit=n_unit, file=fname, status='replace')
        write(n_unit,'(A18, I15)') 'cnt1 =', cnt1
        write(n_unit,'(A18, I15)') 'cnt2 =', cnt2
        write(n_unit,'(A18, I15)') 'cnt3 =', cnt3
        write(n_unit,'(A18, I15)') 'cnt4 =', cnt4
        write(n_unit,'(A18, I15)') 'cnt5 =', cnt5
        write(n_unit,'(A18, I15)') 'cnt6 =', cnt6
        write(n_unit,'(A18, I15)') 'sum =', sum
        do f = 1, cnt1
            write(n_unit,'(A18, F10.4)') 'left(f) =', left(f)
        end do
        do g = 1, cnt2
            write(n_unit,'(A18, F10.4)') 'right(g) =', right(g)
        end do
        do h = 1, cnt3
            write(n_unit,'(A18, F10.4)') 'w(h) =', w(h)
        end do
        do i = 1, cnt4
            write(n_unit,'(A18, F10.4)') 'x(i) =', x(i)
        end do
        do j = 1, cnt5
            write(n_unit,'(A18, F10.4)') 'y(j) =', y(j)
        end do
        do k = 1, cnt6
            write(n_unit,'(A18, F10.4)') 'z(k) =', z(k)
        end do
    
    close(n_unit)

    end if

end subroutine dropletManagement