!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use drop_motion_mod
    implicit none

    integer, save :: phase = 0

    if (dimensional_time(n_time) <= 4.0) then
        if (phase == 0) then
            droplets(250:)%status = -100
            phase = 1
        end if

    else if (dimensional_time(n_time) <= 8.0) then
        if(phase == 1) then
            droplets(250:500)%status = 0
            phase = 2
        end if
    
    else if (dimensional_time(n_time) <= 12.0) then
        if(phase == 2) then
            droplets(500:750)%status = 0
            phase = 3
        end if
    
    else if (dimensional_time(n_time) <= 16.0) then
        if(phase == 3) then
            droplets(750:1000)%status = 0
            phase = 4
        end if
    
    else if (dimensional_time(n_time) <= 20.0) then
        if(phase == 4) then
            droplets(1000:1250)%status = 0
            phase = 5
        end if
    
    else if (dimensional_time(n_time) <= 24.0) then
        if(phase == 5) then
            droplets(1250:1500)%status = 0
            phase = 6
        end if
    
    else if (dimensional_time(n_time) <= 28.0) then
        if(phase == 6) then
            droplets(1500:1750)%status = 0
            phase = 7
        end if

    else if (dimensional_time(n_time) <= 32.0) then
        if(phase == 7) then
            droplets(1750:2000)%status = 0
            phase = 8
        end if
    
    else if (dimensional_time(n_time) <= 36.0) then
        if(phase == 8) then
            droplets(2000:2250)%status = 0
            phase = 9
        end if
    
    else if (dimensional_time(n_time) <= 40.0) then
        if(phase == 9) then
            droplets(2250:2500)%status = 0
            phase = 10
        end if

    else if (dimensional_time(n_time) <= 44.0) then
        if(phase == 10) then
            droplets(2500:2750)%status = 0
            phase = 11
        end if

    else if (dimensional_time(n_time) <= 48.0) then
        if(phase == 11) then
            droplets(2750:3000)%status = 0
            phase = 12
        end if
    
    else if (dimensional_time(n_time) <= 52.0) then
        if(phase == 12) then
            droplets(3000:3250)%status = 0
            phase = 13
        end if
    
    else if (dimensional_time(n_time) <= 56.0) then
        if(phase == 13) then
            droplets(3250:3500)%status = 0
            phase = 14
        end if

    else if (dimensional_time(n_time) <= 60.0) then
        if(phase == 14) then
            droplets(3500:3750)%status = 0
            phase = 15
        end if

    else if (dimensional_time(n_time) <= 64.0) then
        if(phase == 15) then
            droplets(3750:4000)%status = 0
            phase = 16
        end if
    
    else if (dimensional_time(n_time) <= 68.0) then
        if(phase == 16) then
            droplets(4000:4250)%status = 0
            phase = 17
        end if
    
    else if (dimensional_time(n_time) <= 72.0) then
        if(phase == 17) then
            droplets(4250:4500)%status = 0
            phase = 18
        end if

    else if (dimensional_time(n_time) <= 76.0) then
        if(phase == 18) then
            droplets(4500:4750)%status = 0
            phase = 19
        end if

    else if (dimensional_time(n_time) <= 80.0) then
        if(phase == 19) then
            droplets(4750:5000)%status = 0
            phase = 20
        end if
    
    else if (dimensional_time(n_time) <= 84.0) then
        if(phase == 20) then
            droplets(5000:5250)%status = 0
            phase = 21
        end if

    else if (dimensional_time(n_time) <= 88.0) then
        if(phase == 21) then
            droplets(5250:5500)%status = 0
            phase = 22
        end if

    else if (dimensional_time(n_time) <= 92.0) then
        if(phase == 22) then
            droplets(5500:5750)%status = 0
            phase = 23
        end if

    else if (dimensional_time(n_time) <= 96.0) then
        if(phase == 23) then
            droplets(5750:6000)%status = 0
            phase = 24
        end if
    
    else if (dimensional_time(n_time) <= 100.0) then
        if(phase == 24) then
            droplets(6000:6250)%status = 0
            phase = 25
        end if
    
    else if (dimensional_time(n_time) <= 104.0) then
        if(phase == 25) then
            droplets(6250:6500)%status = 0
            phase = 26
        end if

    else if (dimensional_time(n_time) <= 108.0) then
        if(phase == 26) then
            droplets(6500:6750)%status = 0
            phase = 27
        end if

    else if (dimensional_time(n_time) <= 112.0) then
        if(phase == 27) then
            droplets(6750:7000)%status = 0
            phase = 28
        end if
    
    else if (dimensional_time(n_time) <= 116.0) then
        if(phase == 28) then
            droplets(7000:7250)%status = 0
            phase = 29
        end if
    
    else if (dimensional_time(n_time) <= 120.0) then
        if(phase == 29) then
            droplets(7250:7500)%status = 0
            phase = 30
        end if

    else if (dimensional_time(n_time) <= 124.0) then
        if(phase == 30) then
            droplets(7500:7750)%status = 0
            phase = 31
        end if

    else if (dimensional_time(n_time) <= 128.0) then
        if(phase == 31) then
            droplets(7750:8000)%status = 0
            phase = 32
        end if
    
    else if (dimensional_time(n_time) <= 132.0) then
        if(phase == 32) then
            droplets(8000:8250)%status = 0
            phase = 33
        end if
    
    else if (dimensional_time(n_time) <= 136.0) then
        if(phase == 33) then
            droplets(8250:8500)%status = 0
            phase = 34
        end if

    else if (dimensional_time(n_time) <= 140.0) then
        if(phase == 34) then
            droplets(8500:8750)%status = 0
            phase = 35
        end if

    else if (dimensional_time(n_time) <= 144.0) then
        if(phase == 35) then
            droplets(8750:9000)%status = 0
            phase = 36
        end if
    
    else if (dimensional_time(n_time) <= 148.0) then
        if(phase == 36) then
            droplets(9000:9250)%status = 0
            phase = 37
        end if
    
    else if (dimensional_time(n_time) <= 152.0) then
        if(phase == 37) then
            droplets(8250:8500)%status = 0
            phase = 38
        end if

    else if (dimensional_time(n_time) <= 156.0) then
        if(phase == 38) then
            droplets(8500:8750)%status = 0
            phase = 39
        end if

    else if (dimensional_time(n_time) <= 160.0) then
        if(phase == 39) then
            droplets(8750:9000)%status = 0
            phase = 40
        end if

    end if

end subroutine dropletManagement