program dropletCount
    use dropletMotionSimulation
    use boxCounter_m
    use caseName_m
    implicit none
    integer n, num_drop, i_box, num_box, n_unit
    character(255) caseName, fname
    character(:), allocatable :: outFName
    integer, allocatable :: id_array(:), boxCountArray(:)
    type(dropletGroup) dGroup
    type(boxCounter), allocatable :: box_array(:)

    type boxResult_t
        integer num_droplet
        real volume, RoI
    end type
    type(boxResult_t), allocatable :: bResult(:)

    integer, pointer :: nc => nowCase
    integer nc_max

    call case_check(num_case=nc_max)

    DO nc = 1, nc_max                         !実行数だけループ（通常1回）

        caseName = get_caseName()

        call read_and_set_condition(trim(caseName), num_droplet=num_drop)

        box_array = get_box_array(trim(caseName), num_drop)

        num_box = size(box_array)

        if(allocated(boxCountArray)) deallocate(boxCountArray)
        allocate(boxCountArray(num_box))

        fname = trim(caseName)//'/boxTimeSeries.csv'
        open(newunit=n_unit, file=fname, status='replace')

        do n = 0, n_end, outputInterval
            if(n==0) then
                fname = trim(caseName)//'/backup/InitialDistribution.bu'
            else
                write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
            end if

            mainDroplet = read_backup(fname)

            do i_box = 1, num_box
                id_array = mainDroplet%IDinBox(dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
                call box_array(i_box)%add_dropletFlag(id_array)
                boxCountArray(i_box) = size(box_array(i_box)%get_id_array())
            end do

            write(n_unit,'(*(g0:,","))') TimeOnSimu(n), boxCountArray

        end do

        block
            integer i
            real, allocatable :: iniRadDis(:,:), outArray(:,:)

            outFName = trim(caseName)//'/boxInitialRadius.csv'

            if(allocated(bResult)) deallocate(bResult)
            allocate(bResult(num_box))
            do i_box = 1, num_box

                id_array = box_array(i_box)%get_id_array()
                dGroup%droplet = mainDroplet%droplet(id_array)
                bResult(i_box)%num_droplet = size(dGroup%droplet)
                bResult(i_box)%volume = real(dGroup%totalVolume(dim='ml'))
                iniRadDis = dGroup%initialRadiusDistribution()
                if(.not.allocated(outArray)) allocate(outArray(1+num_box, size(iniRadDis, dim=2)))
                outArray(1, :) = iniRadDis(1, :)
                outArray(1+i_box, :) = iniRadDis(2, :)
            end do
            open(newunit=n_unit, file=outFName, status='replace')
            write(n_unit,'(A)') 'radius,B,C,D'
            do i = 1, size(outArray,dim=2)
                write(n_unit,'(*(g0:,","))') outArray(:,i)
            end do
            close(n_unit)
        end block

        bResult(:)%RoI = RateOfInfection(bResult(:)%volume)

        call output_countCSV
        call output_boxVTK

    END DO

    contains

    subroutine output_countCSV
        integer n_unit, i
        character(:), allocatable :: csvFName

        csvFName = trim(caseName)//'/boxCount.csv'
        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_drop,volume[ml],RoI")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_droplet, bResult(i)%volume, bResult(i)%RoI
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

        call mesh%output(trim(caseName)//'/box.vtk', cellScalar=bResult(:)%RoI, scalarName='RoI')

    end subroutine

    elemental real function RateOfInfection(volume)
        real, intent(in) :: volume

        RateOfInfection = 1. - exp(-volume*1.e7 / 900.)

    end function

end program dropletCount