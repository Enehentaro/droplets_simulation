subroutine management_droplets

    use drop_motion_mod

    implicit none

    integer, save :: phase = 0

    if (dimensional_time(n_time) <= 4.0) then

        if (phase == 0) then

        droplets(2000:)%status=-100
        phase = 1

        end if
    
    else if (dimensional_time(n_time) <= 8.0) then

        if (phase == 1) then

        droplets(2000:4000)%status=0
        phase = 2

        end if

    else if (dimensional_time(n_time) <= 12.0) then

        if (phase == 2) then

        droplets(4000:6000)%status=0
        phase = 3

        end if

    else if (dimensional_time(n_time) <= 16.0) then

        if (phase == 3) then

        droplets(6000:8000)%status=0
        phase = 4

        end if

    else if (dimensional_time(n_time) <= 20.0) then

        if (phase == 4) then

        droplets(8000:10000)%status=0
        phase = 5

        end if

    end if

end subroutine management_droplets