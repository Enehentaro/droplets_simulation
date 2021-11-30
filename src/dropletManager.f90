!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use drop_motion_mod
    implicit none

    integer, save :: i = 0
    integer cnt1,cnt2,cnt3,cnt4,sum
    integer j

    if (i == 0) droplets(:)%status = -100

    if (i > 999) return

    do while (Time_onSimulation(n_time) > dble(i)*0.1)
        droplets(10*(i)+1:10*(i+1))%status = 0
        i = i+1    
    end do

    cnt1 = 0
    cnt2 = 0
    cnt3 = 0
    cnt4 = 0
    sum  = 0

    do j = 1, 10000

        if ((-0.036 <= droplets(j)%position(1)) .and. (droplets(j)%position(1) <= 0.015))  then
        if ((-0.468 <= droplets(j)%position(2)) .and. (droplets(j)%position(2) <= -0.447)) then
        if ((0.952  <= droplets(j)%position(3)) .and. (droplets(j)%position(3) <= 0.982))  then
            if (droplets(j)%status == 1) then
                cnt1 = cnt1 + 1
            end if
        end if
        end if
        end if

        if ((-0.020 <= droplets(j)%position(1)) .and. (droplets(j)%position(1) <= 0.031))  then
        if ((-0.468 <= droplets(j)%position(2)) .and. (droplets(j)%position(2) <= -0.433)) then
        if ((0.982  <= droplets(j)%position(3)) .and. (droplets(j)%position(3) <= 1.104))  then
            if (droplets(j)%status == 1) then
                cnt2 = cnt2 + 1
            end if
        end if
        end if
        end if

        if ((-0.014 <= droplets(j)%position(1)) .and. (droplets(j)%position(1) <= 0.014))   then
        if ((-0.468 <= droplets(j)%position(2)) .and. (droplets(j)%position(2) <= -0.3754)) then
        if ((1.104  <= droplets(j)%position(3)) .and. (droplets(j)%position(3) <= 1.1287))  then
            if (droplets(j)%status == 1) then
                cnt3 = cnt3 + 1
            end if
        end if
        end if
        end if

        if ((-0.011 <= droplets(j)%position(1)) .and. (droplets(j)%position(1) <= 0.011))   then
        if ((-0.468 <= droplets(j)%position(2)) .and. (droplets(j)%position(2) <= -0.3601)) then
        if ((1.1287  <= droplets(j)%position(3)) .and. (droplets(j)%position(3) <= 1.163))  then
            if (droplets(j)%status == 1) then
                cnt4 = cnt4 + 1
            end if
        end if
        end if
        end if

    end do

    sum = cnt1 + cnt2 + cnt3 + cnt4

    if ((mod(n_time,interval) == 0)) then
        print*, 'cnt1 =', cnt1
        print*, 'cnt2 =', cnt2
        print*, 'cnt3 =', cnt3
        print*, 'cnt4 =', cnt4
        print*, 'sum  =', sum
    end if    

end subroutine dropletManagement