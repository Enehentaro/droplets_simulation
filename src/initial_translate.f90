program translate
    !!author: Hikaru Konishi
    !!全飛沫に対して、基準点をbefore_dGroupCenter、回転軸をrotation_axisとして、反時計回りにphi[rad]だけ回転させる。
    !!回転後の飛沫中心をafter_dGroupCenterに平行移動させる。
    use virusDroplet_m
    use caseName_m
    implicit none
    type(virusDroplet_t), allocatable :: droplets(:)

    namelist /initial_translate_setting/ fnameDecoration, before_dGroupCenter, after_dGroupCenter, &
    rotation_axis, rotation_angle_deg

    character(50), allocatable :: caseName_array(:)
    character(:), allocatable :: caseName
    character(15) fnameDecoration
    double precision before_dGroupCenter(3), after_dGroupCenter(3)
    character(1) rotation_axis
    double precision rotation_angle_deg

    integer i, n_unit
    double precision, parameter :: PI = acos(-1.d0)
    double precision vec(3), center_displacement(3)
    double precision phi

    ! これから移動を行いたいInitialDistribution.buのあるケース名の取得
    call case_check(caseName_array)
    caseName = trim(caseName_array(1))

    ! optionディレクトリのinitial_translate_setting.nmlの読み込み
    open(newunit = n_unit, file = "option"//"/initial_translate_setting.nml")
        read(n_unit, nml=initial_translate_setting)
    close(n_unit)

    ! 読み込んだ値の確認
    print'(*(g0:))', "fnameDecoration = ", fnameDecoration
    print'(*(g0:," "))', "before_dGroupCenter = ", before_dGroupCenter
    print'(*(g0:," "))', "after_dGroupCenter = ", after_dGroupCenter
    print'(*(g0:))', "rotation_axis = ", rotation_axis
    print'(*(g0:))', "rotation_angle_deg = ", rotation_angle_deg

    ! backupファイルから飛沫の情報を取得
    droplets = read_backup(caseName//"/backup/InitialDistribution.bu")

    ! これから計算を行う上で必要な値の算出
    phi = PI * (rotation_angle_deg / 180.d0) ! [°]から[rad]への変換
    center_displacement = after_dGroupCenter - before_dGroupCenter ! 変位量

    ! 回転軸の切り替え
    select case(rotation_axis)
        case("x")
            ! 回転を行うループ
            do i = 1, size(droplets)
                vec(2) = droplets(i)%position(2) - before_dGroupCenter(2)
                vec(3) = droplets(i)%position(3) - before_dGroupCenter(3)
                droplets(i)%position(2) = cos(phi)*vec(2) - sin(phi)*vec(3) + before_dGroupCenter(2)
                droplets(i)%position(3) = sin(phi)*vec(2) + cos(phi)*vec(3) + before_dGroupCenter(3)
            end do
        case("y")
            do i = 1, size(droplets)
                vec(3) = droplets(i)%position(3) - before_dGroupCenter(3)
                vec(1) = droplets(i)%position(1) - before_dGroupCenter(1)
                droplets(i)%position(3) = cos(phi)*vec(3) - sin(phi)*vec(1) + before_dGroupCenter(3)
                droplets(i)%position(1) = sin(phi)*vec(3) + cos(phi)*vec(1) + before_dGroupCenter(1)
            end do
        case("z")
            do i = 1, size(droplets)      
                vec(1) = droplets(i)%position(1) - before_dGroupCenter(1)
                vec(2) = droplets(i)%position(2) - before_dGroupCenter(2)
                droplets(i)%position(1) = cos(phi)*vec(1) - sin(phi)*vec(2) + before_dGroupCenter(1)
                droplets(i)%position(2) = sin(phi)*vec(1) + cos(phi)*vec(2) + before_dGroupCenter(2)
            end do
    end select

    ! 平行移動を行うループ
    do i = 1, size(droplets)      
        droplets(i)%position(1) = droplets(i)%position(1) + center_displacement(1)
        droplets(i)%position(2) = droplets(i)%position(2) + center_displacement(2)
        droplets(i)%position(3) = droplets(i)%position(3) + center_displacement(3)
    end do

    ! ファイルの出力
    call output_backup(droplets, 'InitialDistribution_'//trim(fnameDecoration)//'.bu')
    call output_droplet_VTK(droplets, 'InitialDistribution_'//trim(fnameDecoration)//'.vtk')
    
end program translate
