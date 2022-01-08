program droplet2CSV
    use dropletGroup_m
    use conditionValue_m
    use dropletEquation_m
    use caseName_m
    implicit none
    integer n, stepEnd, stepInterval
    character(255) caseName, fname
    double precision time
    type(conditionValue_t) condVal
    type(dropletGroup) dGroup

    integer, pointer :: nc => nowCase
    integer nc_max

    call case_check(num_case=nc_max)

DO nc = 1, nc_max                         !実行数だけループ（通常1回）

    caseName = get_caseName()

    call condVal%read(trim(caseName))
    call set_basicVariables_dropletEquation(condVal%dt, condVal%L, condVal%U)

    statusCSV = [0, 1,-1,-2]

    do n = 0, condVal%stepEnd, condVal%outputInterval
        if(n==0) then
            fname = trim(caseName)//'/backup/InitialDistribution.bu'
        else
            write(fname,'("'//trim(caseName)//'/backup/backup_", i0 , ".bu")') n
        end if
        dGroup = read_backup(trim(fname))

        time = TimeOnSimu(step=n, dimension=.true.)
        if(n==0) then
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.true.)
        else
            call dGroup%output_CSV(trim(caseName)//'/particle.csv', time, initial=.false.)
        end if

    end do

END DO
    

end program droplet2CSV