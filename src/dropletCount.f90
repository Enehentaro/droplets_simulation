program dropletCount
    use virusDroplet_m
    use conditionValue_m
    ! use dropletEquation_m
    use boxCounter_m
    use caseName_m
    implicit none
    integer n, i_box, num_box, nc_max
    integer, pointer :: nc => nowCase
    character(255) caseName, fname
    integer, allocatable :: id_array(:)
    type(DropletGroup) mainDroplet, dGroup
    type(conditionValue_t) condVal
    ! type(BasicParameter) baseParam
    type(boxCounter), allocatable :: box_array(:)

    type boxResult_t
        integer num_droplet
        real volume, RoI
    end type
    type(boxResult_t), allocatable :: bResult(:)

    call case_check(num_case = nc_max)
    !print*, 'caseName = ?'
    !read(5, *) caseName

    do nc = 1, nc_max
        caseName = get_caseName()
        call condVal%read(trim(caseName))
        ! baseParam = BasicParameter_(condVal%dt, condVal%L, condVal%U)
    
        box_array = get_box_array(trim(caseName), condVal%num_drop)
    
        num_box = size(box_array)
    
        do n = 0, condVal%stepEnd, condVal%outputInterval
            if(n==0) then
                fname = trim(caseName)//'/backup/InitialDistribution.bu'
            else
                write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
            end if
    
            mainDroplet = read_backup(fname)
    
            do i_box = 1, num_box
                id_array = mainDroplet%IDinBox(dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
                call box_array(i_box)%add_Flag(id_array)
            end do
    
        end do
    
        allocate(bResult(num_box))

        do i_box = 1, num_box
            id_array = box_array(i_box)%get_FlagID()
            dGroup%droplet = mainDroplet%droplet(id_array)
            bResult(i_box)%num_droplet = size(dGroup%droplet)
            bResult(i_box)%volume = real(dGroup%totalVolume() *condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算
        end do

        bResult(:)%RoI = RateOfInfection(bResult(:)%volume) !1分間あたりの感染確率を計算
    
        call output_countCSV
        call output_boxVTK

        deallocate(bResult)

    end do

    contains

    subroutine output_countCSV
        integer n_unit, i
        character(:), allocatable :: csvFName

        csvFName = trim(caseName)//'/BoxCount.csv'
        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_drop,volume[ml],RoI")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_droplet, bResult(i)%volume, bResult(i)%RoI
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK
        use VTK_operator_m
        type(UnstructuredGrid_inVTK) mesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
                                        
        real, allocatable :: xyz(:,:)
        integer, allocatable :: vertices(:,:), types(:)
                                        
        allocate(xyz(3, num_box*8))
        allocate(vertices(8, num_box), types(num_box))
        do i = 1, num_box

            do j = 1, 8
                k = j + 8*(i-1)
                xyz(:,k) = box_array(i)%min_cdn(:) + box_array(i)%width(:)*trans(:,j)
                vertices(j,i) = k
            end do

            types(i) = 11

        end do

        mesh = UnstructuredGrid_inVTK_(xyz, vertices, types)

        call mesh%output(trim(caseName)//'/Box.vtk', cellScalar=bResult(:)%RoI, scalarName='RoI')

    end subroutine

    elemental real function RateOfInfection(volume)
        !1分間あたりの感染確率を計算（もとの資料では1時間あたりの感染確率だが、1分間あたりに換算）
        real, intent(in) :: volume

        RateOfInfection = 1. - exp(-volume*1.e7 / (900./60.))

    end function

end program dropletCount