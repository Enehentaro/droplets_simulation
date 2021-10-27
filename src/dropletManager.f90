!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
! delta_timeが細かいときは,condition.txtの飛沫計算の時間間隔に注意
subroutine management_droplet
    use drop_motion_mod
    implicit none

    integer, parameter          :: num_drop   = 5000
    double precision, parameter :: cough_time = 0.3d0
    double precision, parameter :: delta_time = cough_time / dble(num_drop)
    integer, save               :: cnt        = 1

    if(n_time == 1) then
        droplets(:)%status = -100
    end if

    if(cnt > num_drop) return

    if(dimensional_time(n_time) <= dble(cnt)*delta_time) then
        droplets(cnt)%status = 0
        droplets(cnt + num_drop)%status = 0
        cnt = cnt +1
    end if

end subroutine management_droplet