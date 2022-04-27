program translate
    use dropletGroup_m
    implicit none
    integer n, num_drop, n_unit, i
    character(:), allocatable :: fname, outFName
    type(dropletGroup) dGroup

    dGroup = read_backup('InitialDistribution_convers_A_A2B.bu')

    do i = 1, size(dGroup%droplet)
        ! dGroup%droplet(i)%position(2) = dGroup%droplet(i)%position(2) + 1.36d0  !B_taimen
        ! dGroup%droplet(i)%position(2) = dGroup%droplet(i)%position(2) + 0.3d0  !B2A
        ! dGroup%droplet(i)%position(2) = dGroup%droplet(i)%position(2) - 0.15d0  !fromAtoC
        dGroup%droplet(i)%position(2) = dGroup%droplet(i)%position(2) - 0.075d0  !fromAtoD
    end do

    call dGroup%output_backup(dir='test-test', initial=.true.)
    call dGroup%output_VTK('test-test/test.vtk')
    
end program translate