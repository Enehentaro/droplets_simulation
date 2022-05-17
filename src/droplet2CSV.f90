program droplet2CSV
    use virusDroplet_m
    use conditionValue_m
    use dropletEquation_m
    implicit none
    integer n, stepEnd, stepInterval
    character(255) caseName, fname
    double precision time
    type(conditionValue_t) condVal
    type(BasicParameter) baseParam
    type(DropletGroup) dGroup

    print*, 'caseName = ?'
    read(5, *) caseName

    call condVal%read(caseName)
    baseParam = BasicParameter_(condVal%dt, condVal%L, condVal%U)

    print*, 'End = ?'
    read(5,*) stepEnd

    print*, 'interval = ?'
    read(5,*) stepInterval

    dGroup%statusCSV = [0, 1,-1,-2]

    do n = 0, stepEnd, stepInterval
        write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        dGroup = read_backup(fname)

        time = baseParam%TimeStep2RealTime(step=n, dimension=.true.)
        if(n==0) then
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.true.)
        else
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.false.)
        end if

    end do
    

end program droplet2CSV