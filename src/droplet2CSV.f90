program droplet2CSV
    use dropletMotionSimulation
    implicit none
    integer n, stepEnd, stepInterval
    character(99) caseName, fname
    double precision time

    print*, 'caseName = ?'
    read(5, *) caseName

    call read_and_set_condition(trim(caseName))

    print*, 'End = ?'
    read(5,*) stepEnd

    print*, 'interval = ?'
    read(5,*) stepInterval

    statusCSV = [0, 1,-1,-2]

    do n = 0, stepEnd, stepInterval
        write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        mainDroplet = read_backup(fname)

        time = TimeOnSimu(step=n, dimension=.true.)
        if(n==0) then
            call mainDroplet%output_CSV(trim(caseName)//'/particle.csv', time, initial=.true.)
        else
            call mainDroplet%output_CSV(trim(caseName)//'/particle.csv', time, initial=.false.)
        end if

    end do
    

end program droplet2CSV