program dropletCount
    use dropletMotionSimulation
    use boxCounter_m
    implicit none
    integer n, num_drop, i
    character(99) caseName, fname
    type(boxCounter), allocatable :: box_array(:)
    integer, allocatable :: id_array(:)
    type(dropletGroup) dGroup

    print*, 'caseName = ?'
    read(5, *) caseName

    call read_and_set_condition(trim(caseName), num_droplet=num_drop)

    box_array = get_box_array(trim(caseName), num_drop)

    do n = 0, n_end, outputInterval
        if(n==0) then
            fname = trim(caseName)//'/backup/InitialDistribution.bu'
        else
            write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        end if

        mainDroplet = read_backup(fname)

        do i = 1, size(box_array)
            id_array = mainDroplet%IDinBox(box_array(i)%min_cdn, box_array(i)%max_cdn)
            call box_array(i)%add_dropletFlag(id_array)
        end do

    end do

    call output_countCSV

    contains

    subroutine output_countCSV
        integer n_unit
        character(:), allocatable :: csvFName
        real xyz(3), volume

        csvFName = trim(caseName)//'/boxCount.csv'
        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_drop,volume[ml]")')
            
            do i = 1, size(box_array)
                id_array = box_array(i)%get_id_array()
                dGroup%droplet = mainDroplet%droplet(id_array)
                xyz(:) = real(box_array(i)%min_cdn(:) + box_array(i)%max_cdn(:))/2.0
                volume = real(dGroup%totalVolume(dim='ml'))
                write(n_unit,'(*(g0:,","))') xyz(:), size(dGroup%droplet), volume
            end do

        close(n_unit)

    end subroutine

end program dropletCount