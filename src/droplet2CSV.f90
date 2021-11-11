program droplet2CSV
    use drop_motion_mod
    implicit none
    integer n, stepEnd, stepInterval
    character case_name*99, fname*99

    print*, 'case_name = ?'
    read(5,'(A)') case_name

    print*, 'End = ?'
    read(5,*) stepEnd

    print*, 'interval = ?'
    read(5,*) stepInterval

    statusCSV = [0, 1,-1,-2]

    call first_setting(trim(case_name))

    do n = 0, stepEnd, stepInterval
        write(fname,'("'//trim(case_name)//'/backup/backup", i8.8 , ".bu")') n
        droplets = read_backup(fname)
        call output_droplet_CSV(trim(case_name)//'/particle.csv', droplets(:)%virusDroplet_t, n)

    end do
    

end program droplet2CSV