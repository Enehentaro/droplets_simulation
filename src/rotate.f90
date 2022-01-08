program rotate
    use dropletGroup_m
    implicit none
    integer n, num_drop, n_unit, i
    character(:), allocatable :: fname, outFName
    type(dropletGroup) dGroup
    double precision r, theta, vec(2)
    double precision, parameter :: center(2) = [1.9d0, 1.1d0], PI = acos(-1.d0), &  !A
        phi = atan((1.203d0 - center(2))/(1.838d0 - center(1))) + PI*0.5d0
    ! double precision, parameter :: center(2) = [1.9d0, 0.8d0], PI = acos(-1.d0), &  !B
    !     phi = atan((0.917084d0 - center(2))/(1.85622d0 - center(1))) + PI*0.5d0

    dGroup = read_backup('InitialDistribution_convers_A_A2B.bu')

    print*, 'phi:', phi
    print*, 'PI:', PI

    do i = 1, size(dGroup%droplet)
        vec = dGroup%droplet(i)%position(1:2) - center
        r = norm2(vec)
        theta = atan(vec(2)/vec(1))
        if (theta < 0.d0) theta = theta + PI
        ! print*, 'theta:', theta
        dGroup%droplet(i)%position(1) = r*cos(theta+phi) + center(1)
        dGroup%droplet(i)%position(2) = r*sin(theta+phi) + center(2)
    end do

    call dGroup%output_backup(dir='test-test', initial=.true.)
    call dGroup%output_VTK('test-test/test.vtk')
    
end program rotate