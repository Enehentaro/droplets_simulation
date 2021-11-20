!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use drop_motion_mod
    implicit none

    integer, save :: i = 0

    if (i == 0) droplets(:)%status = -100

    if (i > 999) return

    do while (Time_onSimulation(n_time) > dble(i)*0.1)
        droplets(10*(i)+1:10*(i+1))%status = 0
        i = i+1    
    end do

end subroutine dropletManagement