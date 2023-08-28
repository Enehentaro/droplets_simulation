program dropletCount
    !!ボックスを通過した飛沫をカウントする。
    !!飛沫計算の出力ファイルを順に読み込み、各ボックスに対して内外判定を行う。
    !!ボックス側では、通過した飛沫のIDしか見ておらず、同じIDの飛沫のダブルカウントなどは起こらない。
    use virusDroplet_m
    use conditionValue_m
    ! use dropletEquation_m
    use boxCounter_m
    use caseName_m
    implicit none
    integer n, i_box, num_box, caseID , i 
    character(50), allocatable :: caseName_array(:)
    character(:), allocatable :: caseName, fname
    integer, allocatable :: id_array(:)
    type(virusDroplet_t), allocatable :: mainDroplets(:), droplets(:)
    type(conditionValue_t) condVal
    ! type(BasicParameter) baseParam
    type(boxCounter), allocatable :: box_array(:)

    type boxResult_t
        integer num_droplet
        real volume, RoI
        integer num_adhes
    end type
    type(boxResult_t), allocatable :: bResult(:)

    call case_check(caseName_array)
    !print*, 'caseName = ?'
    !read(5, *) caseName

    do caseID = 1, size(caseName_array)                     
            !ケースの数だけ繰り返す。複数の場合はテキストファイルに出力しておく必要がある
        caseName = trim(caseName_array(caseID))             
            !caseName_arrayにある名前には余分な空白があるからtrimで文字の部分だけを取り出している　
        condVal = read_condition(caseName)
        ! baseParam = BasicParameter_(condVal%dt, condVal%L, condVal%U)
    
        box_array = get_box_array(caseName, condVal%num_drop)
    
        num_box = size(box_array)   
            !データのサイズだから今であればboxの数
    
        ! do n = 0, condVal%stepEnd, condVal%outputInterval
        !     if(n==0) then
        !         fname = caseName//'/backup/InitialDistribution.bu'
        !     else
        !         block
        !             character(255) str
        !             write(str,'("'//caseName//'/backup/backup_", i0 , ".bu")') n
        !             fname = trim(str)
        !         end block
        !     end if
    
        !     mainDroplets = read_backup(fname)
    
        !     do i_box = 1, num_box
        !         id_array = dropletIDinBox(mainDroplets, dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
        !         call box_array(i_box)%add_Flag(id_array)
        !     end do
    
        ! end do

        fname = caseName//'/backup/backup_10000.bu'
        mainDroplets = read_backup(fname)
            do i_box = 1, num_box
                id_array = dropletIDinBox(mainDroplets, dble(box_array(i_box)%min_cdn), dble(box_array(i_box)%max_cdn))
                call box_array(i_box)%add_Flag(id_array)
            end do
            !浮遊している飛沫にフラグを立てる

        allocate(bResult(num_box))

        do i_box = 1, num_box
            id_array = box_array(i_box)%get_FlagID()
            droplets = mainDroplets(id_array)
            bResult(i_box)%num_droplet = size(droplets)
            bResult(i_box)%volume = real(dropletTotalVolume(droplets) *condVal%L**3 * 1.d6 )    !有次元化[m^3]したのち、[ml]に換算
            bResult(i_box)%num_adhes = 0
            do i = 1 , size(droplets)
                if(.not.droplets(i)%isFloating()) then
                    bResult(i_box)%num_adhes = bResult(i_box)%num_adhes + 1 
                endif
            enddo
            !adhes は付着した飛沫の数であるから、その総数を計算している
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

        csvFName = caseName//'/BoxCount.csv'
        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,num_drop,volume[ml],RoI,num_adhes")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%num_droplet,&
                bResult(i)%volume, bResult(i)%RoI,bResult(i)%num_adhes
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK
        use VTK_operator_m
        type(UnstructuredGrid_inVTK) mesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
                                        
        real, allocatable :: xyz(:,:)
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

        call mesh%output(caseName//'/Box.vtk', cellScalar=bResult(:)%RoI, scalarName='RoI')

    end subroutine

    !>1分間あたりの感染確率を計算（もとの資料では1時間あたりの感染確率だが、1分間あたりに換算）
    elemental real function RateOfInfection(volume)
        real, intent(in) :: volume

        RateOfInfection = 1. - exp(-volume*1.e7 / (900./60.))

    end function

end program dropletCount