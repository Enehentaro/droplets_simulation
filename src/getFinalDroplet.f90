program dropletCount
    !!ボックスを通過した飛沫をカウントする。
    !!飛沫計算の出力ファイルを順に読み込み、各ボックスに対して内外判定を行う。
    !!ボックス側では、通過した飛沫のIDしか見ておらず、同じIDの飛沫のダブルカウントなどは起こらない。
    use virusDroplet_m
    use conditionValue_m
    ! use dropletEquation_m
    use boxCounter_m
    ! use caseName_m
    use simpleFile_reader
    use path_operator_m
    implicit none
    integer caseID, officeID
    integer endStep
    character(255) path2officelist, outputFName, backupFName
    character(50), allocatable :: office_array(:)
    character(50), allocatable :: caseName_array(:)
    character(:), allocatable :: caseName, path2mainDir, officeName, path2caseDir
    integer, allocatable :: id_array(:)
    type(virusDroplet_t), allocatable :: mainDroplets(:), droplets(:)
    type(conditionValue_t) condVal
    ! integer outputInterval
    type(boxCounter), allocatable :: box_array(:)
    real deltaTime

    type result_t
        integer num_floating
        real volume
    end type

    type(result_t), allocatable :: result_overCases(:)

    print*, 'path2officelist = ?'
    read(5, '(A)') path2officelist

    call read_textRecord(trim(path2officelist), office_array)

    call get_DirFromPath(trim(path2officelist), path2mainDir)

    do officeID = 1, size(office_array)

        officeName = trim(office_array(officeID))

        path2caseDir = trim(path2mainDir) // officeName // '_data_droplets/'

        ! call case_check(caseName_array)
        call read_textRecord(path2caseDir//'case_list.txt', caseName_array)
        allocate(result_overCases(size(caseName_array)))

        endStep = 600000
            
        do caseID = 1, size(caseName_array)
            caseName = path2caseDir // trim(caseName_array(caseID))
            condVal = read_condition(caseName)
            deltaTime = real(condVal%dt * condVal%L/condVal%U)
        
            write(backupFName,'("'//caseName//'/backup/backup_", i0 , ".bu")') endStep
            mainDroplets = read_backup(trim(backupFName))

            id_array = dropletIDinState(mainDroplets, "floating")
            droplets = mainDroplets(id_array)
            result_overCases(caseID)%num_floating = size(droplets)
            result_overCases(caseID)%volume = real(dropletTotalVolume(droplets) * condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算

        end do

            
        write(outputFName, '("floating_count_results/", A, ".csv")') officeName
        call output_CSV_overCases(trim(outputFName))

        deallocate(result_overCases)

    end do

    contains

    ! subroutine calcRoI_and_output
    !     use path_operator_m
    !     real erapsedTime
    !         !! 経過時間 [ h ]

    !     type(boxResult_t), allocatable :: bResult(:)

    !     character(255) :: filename, filename2
    !     character(:), allocatable :: output_path

    !     erapsedTime = n*deltaTime / 3600.   ! 経過時間 [ h ]
        
    !     allocate(bResult(num_box))

    !     do i_box = 1, num_box
    !         id_array = result_overCases(caseID)%box%get_FlagID()
    !         dGroup%droplet = mainDroplet%droplet(id_array)
    !         bResult(i_box)%num_droplet = size(dGroup%droplet)
    !         bResult(i_box)%volume = real(dGroup%totalVolume() * condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算
    !     end do

    !     bResult(:)%RoI = RateOfInfection(bResult(:)%volume, erapsedTime) !1分間あたりの感染確率を計算
    
    !     output_path = caseName//'/BoxCount_from1sec'
    !     call make_directory(output_path)
    !     ! call remove_directory(caseName//'/BoxCount')

    !     write(filename, '("'//output_path//'/box_", i0, ".csv")') n
    !     call output_countCSV(trim(filename), bResult)

    !     write(filename, '("'//output_path//'/Box_", i0, ".vtk")') n
    !     write(filename2, '("'//output_path//'/Box_c_", i0, ".vtk")') n
    !     call output_boxVTK(trim(filename), trim(filename2), bResult)

    ! end subroutine

    subroutine output_CSV_overCases(csvFName)
        character(*), intent(in) :: csvFName
        integer n_unit, i!, iost
        ! character(128) errmsg
        ! logical  :: isitopened

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
            write(n_unit, '("casename,num_floating,volume[ml]")')
            
            do i = 1, size(caseName_array)
                write(n_unit,'(*(g0:,","))') trim(caseName_array(i)), result_overCases(i)%num_floating, result_overCases(i)%volume
            end do

        close(n_unit)

    end subroutine

    subroutine output_countCSV(csvFName, bResult)
        character(*), intent(in) :: csvFName
        type(result_t), intent(in) :: bResult(:)
        integer n_unit, i

        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_floating,volume[ml]")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_floating, bResult(i)%volume
            end do

        close(n_unit)

    end subroutine

    ! subroutine output_boxVTK(vtkFName, vtkFName2, bResult)
    !     use VTK_operator_m
    !     character(*), intent(in) :: vtkFName
    !     character(*), intent(in) :: vtkFName2
    !     type(boxResult_t), intent(in) :: bResult(:)
    !     type(UnstructuredGrid_inVTK) mesh
    !     integer i, j, k
    !     real, parameter :: trans(3,8) = reshape([ &
    !                                         0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
    !                                         0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
                                        
    !     real, allocatable :: xyz(:,:), num_droplet(:)
    !     integer, allocatable :: vertices(:,:), types(:)
                                        
    !     allocate(xyz(3, num_box*8))
    !     allocate(vertices(8, num_box), types(num_box))
    !     do i = 1, num_box

    !         do j = 1, 8
    !             k = j + 8*(i-1)
    !             xyz(:,k) = box_array(i)%min_cdn(:) + box_array(i)%width(:)*trans(:,j)
    !             vertices(j,i) = k
    !         end do

    !         types(i) = 11

    !     end do

    !     mesh = UnstructuredGrid_inVTK_(xyz, vertices, types)

    !     call mesh%output(vtkFName, cellScalar=bResult(:)%RoI, scalarName='RoI')

    !     num_droplet = real(bResult(:)%num_droplet)
    !     call mesh%output(vtkFName2, cellScalar=num_droplet, scalarName='Droplets')

    ! end subroutine

    
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