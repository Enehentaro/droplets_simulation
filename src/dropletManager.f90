!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use dropletMotionSimulation
    implicit none
    ! type(dropletGroup) dGroup
    double precision, save :: next_time = 0.d0
    integer, save :: i = 0

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

    if(timeStep==0) mainDroplet%droplet(:)%status = -99

    if(TimeOnSimu(dimension=.true.) >= next_time) then
        ! dGroup = generate_dropletGroup(150)
        ! call mainDroplet%append(dGroup)
        mainDroplet%droplet(16*i + 1 : 16*(i+1))%status = 0
        i = i + 1
        next_time = next_time + 0.1d0
    end if
    
end subroutine dropletManagement