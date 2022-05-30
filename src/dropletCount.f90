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
        ! addition
        integer adherent_droplet
        integer float_droplet
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
    
        ! do n = 0, condVal%stepEnd, condVal%outputInterval
        !     if(n==0) then
        !         fname = trim(caseName)//'/backup/InitialDistribution.bu'
        !     else
        !         write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        !     end if
    
        !     mainDroplet = read_backup(fname)
    
        !     do i_box = 1, num_box
        !         id_array = mainDroplet%IDinBox(dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
        !         call box_array(i_box)%add_Flag(id_array)
        !     end do
    
        ! end do

        ! addition
        fname = trim(caseName)//'/backup/backup_6000000.bu'
        
        mainDroplet = read_backup(fname)
    
            do i_box = 1, num_box
                id_array = mainDroplet%IDinBox(dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
                call box_array(i_box)%add_Flag(id_array)
            end do        
    
        allocate(bResult(num_box))

        do i_box = 1, num_box
            id_array = box_array(i_box)%get_FlagID()
            dGroup%droplet = mainDroplet%droplet(id_array)
            bResult(i_box)%num_droplet = size(dGroup%droplet)
            bResult(i_box)%volume = real(dGroup%totalVolume() *condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算
            ! addition
            bResult(i_box)%adherent_droplet= count(dGroup%droplet(:)%status == 1)
            bResult(i_box)%float_droplet= count(dGroup%droplet(:)%status == 0)
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

            ! addition
            write(n_unit, '("x,y,z,num_drop,volume[ml],RoI,adherent_drop,float_drop")')
            
            do i = 1, size(box_array)
                ! addition
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_droplet, bResult(i)%volume, bResult(i)%RoI,&
                bResult(i)%adherent_droplet, bResult(i)%float_droplet
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK
        use vtkMesh_operator_m
        type(vtkMesh) mesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
                                        
        call mesh%allocation_node(num_box*8)
        call mesh%allocation_cell(num_box)

        do i = 1, num_box

            do j = 1, 8
                k = 8*(i-1) + j - 1
                mesh%node_array(k)%coordinate(:) = box_array(i)%min_cdn(:) + box_array(i)%width(:)*trans(:,j)
            end do

            mesh%cell_array(i-1)%nodeID = [(8*(i-1) + j - 1, j = 1, 8)]
            mesh%cell_array(i-1)%n_TYPE = 11

        end do

        call mesh%output(trim(caseName)//'/Box.vtk', cellScalar=bResult(:)%RoI, scalarName='RoI')

    end subroutine

    elemental real function RateOfInfection(volume)
        !1分間あたりの感染確率を計算（もとの資料では1時間あたりの感染確率だが、1分間あたりに換算）
        real, intent(in) :: volume

        RateOfInfection = 1. - exp(-volume*1.e7 / (900./60.))

    end function

end program dropletCount