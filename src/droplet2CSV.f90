program droplet2CSV
    use dropletGroup_m
    use conditionValue_m
    use dropletEquation_m
    implicit none
    integer n, stepEnd, stepInterval
    character(255) caseName, fname
    double precision time
    type(conditionValue_t) condVal
    type(dropletGroup) dGroup

    print*, 'caseName = ?'
    read(5, *) caseName

    call condVal%read(caseName)
    call set_basicVariables_dropletEquation(condVal%dt, condVal%L, condVal%U)

    print*, 'End = ?'
    read(5,*) stepEnd

    print*, 'interval = ?'
    read(5,*) stepInterval

    statusCSV = [0, 1,-1,-2]

    do n = 0, stepEnd, stepInterval
        write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        dGroup = read_backup(fname)

        time = TimeOnSimu(step=n, dimension=.true.)
        if(n==0) then
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.true.)
        else
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.false.)
        end if

    end do
    

end program droplet2CSV