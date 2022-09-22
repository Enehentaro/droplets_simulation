program dropletCount
    !!ボックスを通過した飛沫をカウントする。
    !!飛沫計算の出力ファイルを順に読み込み、各ボックスに対して内外判定を行う。
    !!ボックス側では、通過した飛沫のIDしか見ておらず、同じIDの飛沫のダブルカウントなどは起こらない。
    use virusDroplet_m
    use conditionValue_m
    ! use dropletEquation_m
    use boxCounter_m
    use caseName_m
    use simpleFile_reader
    use path_operator_m
    implicit none
    integer n, i_box, num_box, caseID, startStep
    character(255) path2caselist, outputFName
    character(50), allocatable :: caseName_array(:)
    character(:), allocatable :: caseName, fname, path2mainDir
    integer, allocatable :: id_array(:)
    type(DropletGroup) mainDroplet, dGroup
    type(conditionValue_t) condVal
    ! integer outputInterval
    type(boxCounter), allocatable :: box_array(:)
    real deltaTime
    double precision min_cdn(3), max_cdn(3)

    type boxResult_t
        integer num_droplet
        real volume, RoI
    end type

    type(boxResult_t), allocatable :: results(:)

    print*, 'path2caselist = ?'
    read(5, '(A)') path2caselist

    print*, 'startStep = ?'
    read(5, *) startStep

    print*, 'outputFName = ?'
    read(5, '(A)') outputFName

    call get_DirFromPath(trim(path2caselist), path2mainDir)
    ! call case_check(caseName_array)
    call read_textRecord(trim(path2caselist), caseName_array)

    allocate(results(size(caseName_array)))

    do caseID = 1, size(caseName_array)
        caseName = path2mainDir // trim(caseName_array(caseID))
        condVal = read_condition(caseName)
        deltaTime = real(condVal%dt * condVal%L/condVal%U)
    
        box_array = get_box_array(path2mainDir // 'box.csv', condVal%num_drop)
    
        num_box = size(box_array)
    
        do n = startStep, condVal%stepEnd, condVal%outputInterval
            if(n==0) then
                fname = caseName//'/backup/InitialDistribution.bu'
            else
                block
                    character(255) str
                    write(str,'("'//caseName//'/backup/backup_", i0 , ".bu")') n
                    fname = trim(str)
                end block
            end if
    
            mainDroplet = read_backup(fname)
    
            do i_box = 1, num_box
                min_cdn = dble(box_array(i_box)%min_cdn)
                max_cdn = dble(box_array(i_box)%max_cdn)
                id_array = mainDroplet%IDinBox(min_cdn, max_cdn)
                call box_array(i_box)%add_Flag(id_array)
            end do

            ! if(mod(n, outputInterval) == 0 .and. n /= 0) call calcRoI_and_output
    
        end do

        id_array = box_array(1)%get_FlagID()
        results(caseID)%num_droplet = size(id_array)
        dGroup%droplet = mainDroplet%droplet(id_array)
        results(caseID)%volume = real(dGroup%totalVolume() * condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算

    end do

    call output_CSV_overCases(outputFName)

    contains

    subroutine calcRoI_and_output
        use path_operator_m
        real erapsedTime
            !! 経過時間 [ h ]

        type(boxResult_t), allocatable :: bResult(:)

        character(255) :: filename, filename2
        character(:), allocatable :: output_path

        erapsedTime = n*deltaTime / 3600.   ! 経過時間 [ h ]
        
        allocate(bResult(num_box))

        do i_box = 1, num_box
            id_array = box_array(i_box)%get_FlagID()
            dGroup%droplet = mainDroplet%droplet(id_array)
            bResult(i_box)%num_droplet = size(dGroup%droplet)
            bResult(i_box)%volume = real(dGroup%totalVolume() * condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算
        end do

        bResult(:)%RoI = RateOfInfection(bResult(:)%volume, erapsedTime) !1分間あたりの感染確率を計算
    
        output_path = caseName//'/BoxCount_from1sec'
        call make_directory(output_path)
        ! call remove_directory(caseName//'/BoxCount')

        write(filename, '("'//output_path//'/box_", i0, ".csv")') n
        call output_countCSV(trim(filename), bResult)

        write(filename, '("'//output_path//'/Box_", i0, ".vtk")') n
        write(filename2, '("'//output_path//'/Box_c_", i0, ".vtk")') n
        call output_boxVTK(trim(filename), trim(filename2), bResult)

    end subroutine

    subroutine output_CSV_overCases(csvFName)
        character(*), intent(in) :: csvFName
        integer n_unit, i, iost
        character(128) errmsg
        logical  :: isitopened

        print*, 'output: ', csvFName

        ! inquire(file=csvFName, exist=isitopened)
        ! if(isitopened) then
        !     print*, "We are there"
        !     call chmod(csvFName, "u+wrx")
        !     open(newunit=n_unit, file=csvFName, status='old', action = "write")
        ! else
        !     print*, "not there"
        !     call chmod(csvFName, "u+wrx")
        !     open(newunit=n_unit, file=csvFName, status='new', action = "write")
        ! end if

        ! open(unit=220, file="~/test2.txt")

        open(newunit=n_unit, file=csvFName, status='replace', action = "write")
        ! open(newunit=n_unit, file=csvFName, status='replace', action = "write", iostat = iost, iomsg = errmsg)
            ! if ( iost /= 0 ) then
            !     print *, "cannot open "//csvFName, iost
            !     error stop "ERROR :: "//trim(errmsg)
            ! end if
            write(n_unit, '("casename,num_drop,volume[ml]")')
            
            do i = 1, size(caseName_array)
                write(n_unit,'(*(g0:,","))') trim(caseName_array(i)), results(i)%num_droplet, results(i)%volume
            end do

        close(n_unit)

    end subroutine

    subroutine output_countCSV(csvFName, bResult)
        character(*), intent(in) :: csvFName
        type(boxResult_t), intent(in) :: bResult(:)
        integer n_unit, i

        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_drop,volume[ml],RoI")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_droplet, bResult(i)%volume, bResult(i)%RoI
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK(vtkFName, vtkFName2, bResult)
        use VTK_operator_m
        character(*), intent(in) :: vtkFName
        character(*), intent(in) :: vtkFName2
        type(boxResult_t), intent(in) :: bResult(:)
        type(UnstructuredGrid_inVTK) mesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
                                        
        real, allocatable :: xyz(:,:), num_droplet(:)
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

        call mesh%output(vtkFName, cellScalar=bResult(:)%RoI, scalarName='RoI')

        num_droplet = real(bResult(:)%num_droplet)
        call mesh%output(vtkFName2, cellScalar=num_droplet, scalarName='Droplets')

    end subroutine

    
    elemental real function RateOfInfection(volume, erapsedTime)
        !! 感染確率を計算（もとの資料では1時間あたりの感染確率だが、経過時間あたりに換算）

        real, intent(in) :: volume
            !! 飛沫総体積 [ ml ]

        real, intent(in) :: erapsedTime
            !! 経過時間 [ h ]

        real, parameter :: n_v = 1.e7
            !! ウイルス密度 [ viral copy / ml ]

        real, parameter :: N_0 = 900.
            !! 感染に至るウイルス量の指標 [ viral copy / h ]

        RateOfInfection = 1. - exp( - volume * n_v / (N_0 * erapsedTime))

    end function

end program dropletCount