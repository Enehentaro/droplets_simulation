!>全飛沫を、center(x,y)中心にphi[rad]だけ回転させ(z方向には回転しない),全飛沫の中心center(x,y,z)をcenter_translated(X,Y,Z)に移動させるサンプルコード
!> by Hikaru Konishi
program translate
    use virusDroplet_m
    implicit none
    character(15) infected_person
    integer i
    type(dropletGroup) dGroup
    double precision vec(2)
    double precision, parameter :: center(3) = [4.73d0, 2.18d0, 1.255d0], PI = acos(-1.d0), phi = PI*0.d0
    double precision center_displacement(3)
    double precision, parameter :: center_translated(3) = [0.955d0, 0.54d0, 1.255d0]

    print*, 'Who is infected ?'
    read(5, '(A)') infected_person

    vec = 0.0d0
    center_displacement = 0.0d0
    center_displacement = center_translated - center

    dGroup = read_backup('InitialDistribution.bu')    !読み込むBUファイル名

    print*, 'phi:', phi
    print*, 'PI:', PI

    do i = 1, size(dGroup%droplet)      !回転を行うループ
        vec = dGroup%droplet(i)%position(1:2) - center(1:2)
        dGroup%droplet(i)%position(1) = cos(phi)*vec(1) - sin(phi)*vec(2) + center(1)
        dGroup%droplet(i)%position(2) = sin(phi)*vec(1) + cos(phi)*vec(2) + center(2)
    end do

    do i = 1, size(dGroup%droplet)      !平行移動を行うループ
        dGroup%droplet(i)%position(1) = dGroup%droplet(i)%position(1) + center_displacement(1)
        dGroup%droplet(i)%position(2) = dGroup%droplet(i)%position(2) + center_displacement(2)
        dGroup%droplet(i)%position(3) = dGroup%droplet(i)%position(3) + center_displacement(3)
    end do

    call dGroup%output_backup('InitialDistribution_'//trim(infected_person)//'.bu')  !BUファイル出力
    call dGroup%output_VTK('InitialDistribution_'//trim(infected_person)//'.vtk')                !確認VTKファイル出力
    
end program translate