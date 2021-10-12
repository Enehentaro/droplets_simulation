module unstructuredGrid_mod
    implicit none
    character(3), private :: FILE_TYPE  !ファイル形式

    integer, allocatable :: ICN(:,:)                !要素所有節点ID
    integer, allocatable :: NoB(:)                  !要素所有境界面の数
    integer, allocatable :: ICB(:,:)                !要素所有境界面ID
    integer, allocatable :: NBN(:,:)                !境界面所有節点ID
    
    integer, allocatable :: CELL_TYPE(:)   !要素タイプ（テトラ0、プリズム1、ピラミッド2）
    integer, allocatable, private :: NUM_NC(:)               !隣接要素数
    integer, allocatable, private :: NEXT_CELL(:,:)          !隣接要素ID

    double precision, allocatable :: CDN(:,:)       !節点座標
    double precision, allocatable :: VELC(:,:)      !要素流速
    double precision, allocatable, private :: CENC(:,:)      !要素重心
    double precision, allocatable, private :: WIDC(:)        !要素の1辺長さ
    double precision, allocatable :: CENF(:,:)    !面重心
    double precision, allocatable :: MOVF(:,:)    !面重心移動量
    double precision, allocatable :: NVECF(:,:)     !面法線ベクトル

    interface read_unstructuredGrid
        module procedure read_unstructuredGrid_byNAME
        module procedure read_unstructuredGrid_byNumber
    end interface

    contains

    subroutine check_FILE_TYPE(FNAME)
        character(*), intent(in) :: FNAME

        if(index(FNAME, '.vtk') > 0) then
            FILE_TYPE = 'VTK'

        else if(index(FNAME, '.inp') > 0) then
            FILE_TYPE = 'INP'

        else if(index(FNAME, '.fld') > 0) then
            FILE_TYPE = 'FLD'

        else
            print*, 'FILE_TYPE NG:', FNAME
            stop
        end if

        print*, 'FILE_TYPE:', FILE_TYPE

    end subroutine check_FILE_TYPE

    subroutine read_unstructuredGrid_byNAME(FNAME)
        character(*), intent(in) :: FNAME
        integer num_cell

        select case(FILE_TYPE)
            case('VTK')
                call read_VTK(FNAME)

            case('INP')
                call read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('FLD')
                call read_FLD(FNAME)

            case default
                print*,'FILE_TYPE NG:', FILE_TYPE
                STOP
                    
        end select

        num_cell = size(ICN(:,:), dim=2)

        if(.not.allocated(CENC)) allocate(CENC(3, num_cell), WIDC(num_cell))
            
    end subroutine read_unstructuredGrid_byNAME

    subroutine read_unstructuredGrid_byNumber(path, digits_fmt, FNUM)
        character(*), intent(in) :: path
        integer, intent(in) :: FNUM
        character(99) :: FNAME
        character(4) digits_fmt
        integer num_cell

        select case(FILE_TYPE)
            case('VTK')

                write(FNAME,'("'//trim(path)//'",'//digits_fmt//',".vtk")') FNUM
                call read_VTK(FNAME)

            case('INP')

                if(FNUM==0) then
                    write(FNAME,'("'//trim(path)//'",'//digits_fmt//',".inp")') 1
                else
                    write(FNAME,'("'//trim(path)//'",'//digits_fmt//',".inp")') FNUM
                end if

                call read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('FLD')

                if (FNUM <= 9) then
                    write(FNAME,'("'//trim(path)//'",i1.1,".fld")') FNUM

                else if (FNUM <= 99) then
                    write(FNAME,'("'//trim(path)//'",i2.2,".fld")') FNUM

                else if(FNUM <= 999) then
                    write(FNAME,'("'//trim(path)//'",i3.3,".fld")') FNUM

                else if(FNUM <= 9999) then
                    write(FNAME,'("'//trim(path)//'",i4.4,".fld")') FNUM

                else
                    write(FNAME,'("'//trim(path)//'",i5.5,".fld")') FNUM

                end if

                call read_FLD(FNAME)

            case default
                print*,'FILE_TYPE NG:', FILE_TYPE
                STOP
                
        end select

        num_cell = size(ICN(:,:), dim=2)

        if(.not.allocated(CENC)) allocate(CENC(3, num_cell), WIDC(num_cell))
            
    end subroutine read_unstructuredGrid_byNumber

    subroutine read_VTK(FNAME, pointdata)
        character(*), intent(in) :: FNAME
        logical, optional :: pointdata
        integer II,KK,IIH, n_unit, KKMX, IIMX
        character AAA*7
        double precision, allocatable :: UVWK(:,:)
        logical prism_flag

        prism_flag = .false.

        print*, 'READ_VTK:', FNAME
            
        open(newunit=n_unit,FILE=FNAME, STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,*) AAA,KKMX
                
            if(.not.allocated(CDN)) allocate(CDN(3,KKMX), source=0.0d0)      
            DO KK = 1, KKMX
                read(n_unit,*) CDN(1,KK),CDN(2,KK),CDN(3,KK)
            END DO
            read(n_unit,'()')
            read(n_unit,'(A,I12)') AAA,IIMX
            
            if(.not.allocated(ICN)) allocate(ICN(6,IIMX), source=0)
            if(.not.allocated(CELL_TYPE)) allocate(CELL_TYPE(IIMX), source=-1)      
            DO II = 1, IIMX
                read(n_unit,fmt='(I12)',advance='no') IIH
                IF(IIH==6) THEN
                    read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II),ICN(5,II),ICN(6,II)
                    CELL_TYPE(II)=  1
                    prism_flag = .true.
                ELSE IF(IIH==5) THEN
                    read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II),ICN(5,II)
                    CELL_TYPE(II) = 2
                ELSE IF(IIH==4) THEN
                    read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II)
                    CELL_TYPE(II) = 0
                END IF
            END DO
            read(n_unit,'()')
            
            
            ICN(:,:) = ICN(:,:) + 1

            if(.not.allocated(VELC)) allocate(VELC(3,IIMX))

            read(n_unit,'()')  !CELL_TYPES
            DO II = 1, IIMX
                read(n_unit,'()')
            END DO
            read(n_unit,'()')

            if(.not.prism_flag) then
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                DO II = 1, IIMX
                    read(n_unit,'()')
                END DO
                read(n_unit,'()')
            
                DO II = 1, IIMX
                    read(n_unit,*) VELC(1,II),VELC(2,II),VELC(3,II)
                END DO

            else

                if(present(pointdata)) then
                    if(pointdata) then
                        allocate(UVWK(3,KKMX))
                        read(n_unit,'()')
                        read(n_unit,'()')
                        DO KK = 1, KKMX
                        read(n_unit,*) UVWK(1,KK),UVWK(2,KK),UVWK(3,KK)
                        END DO
                        call point2cell(UVWK, VELC, ICN, CELL_TYPE)
                        close(n_unit)
                        return
                    end if
                end if
                
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                DO II = 1, IIMX
                    read(n_unit,'()')
                END DO
                read(n_unit,'()')
                read(n_unit,'()')
                DO II = 1, IIMX
                    read(n_unit,*) AAA
                END DO
                read(n_unit,'()')
            
                DO II = 1, IIMX
                    read(n_unit,*) VELC(1,II),VELC(2,II),VELC(3,II)
                END DO
            end if

            
        close(n_unit)
            
    end subroutine read_VTK

    subroutine read_INP(FNAME)
        !  INPファイルを読み込み、節点データを要素データに変換する
        character(*), intent(in) :: FNAME
        INTEGER II,II2,KK,AAmax, IIMX2, AA, n_unit
        integer :: IITETMX, IIPRSMX, IIPYRMX, KKMX, IIMX
        character(6) cellshape
        double precision, allocatable :: UVWK(:,:)
        integer, allocatable :: ICN2(:,:)
        integer, allocatable :: CELL_TYPE2(:)

        print*, 'READ_INP:', trim(FNAME)
            
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
            read(n_unit,*)KKMX,IIMX2
            print*,'KKMX,IIMX2=',KKMX,IIMX2
            
            if(.not.allocated(CDN)) allocate(CDN(3,KKMX), source=0.0d0)
            allocate(ICN2(6,IIMX2), source=0)
            allocate(CELL_TYPE2(IIMX2), source=0)
                
            DO KK = 1, KKMX
                read(n_unit,*)AA,CDN(1,KK),CDN(2,KK),CDN(3,KK)
            END DO
                
            II = 0
            IITETMX = 0
            IIPRSMX = 0
            IIPYRMX = 0
            DO II2 = 1, IIMX2
    
                read(n_unit,fmt='(I10)',advance='no')AA  !ここはセル番号なので無視
                read(n_unit,fmt='(I6)',advance='no')AA  !ここはなんかしらん（だいたいゼロ）
                read(n_unit,fmt='(A6)',advance='no')cellshape
    
                cellshape = adjustl(cellshape)    !左詰め
                IF ((cellshape=='tet').or.(cellshape=='prism').or.(cellshape=='pyr')) THEN
                    II = II +1
        
                    if (cellshape=='tet') then
                        CELL_TYPE2(II) = 0
                        IITETMX = IITETMX +1
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II)    
                        if (ICN2(1,II)==0.or.ICN2(4,II)==0) print*, 'ICN2_WARNING_tet:', ICN2(:,II)
                    
                    ELSE IF(cellshape=='prism') THEN
                        CELL_TYPE2(II) = 1
                        IIPRSMX = IIPRSMX +1
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II),ICN2(5,II),ICN2(6,II)
                        if (ICN2(1,II)==0.or.ICN2(6,II)==0) print*, 'ICN2_WARNING_prism:', ICN2(:,II)
        
                    ELSE IF(cellshape=='pyr') THEN
                        CELL_TYPE2(II) = 2
                        IIPYRMX = IIPYRMX +1
                        read(n_unit,*)ICN2(5,II),ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II) !INPは最初が山頂点であり、VTKでは最後が山頂点のため、読み込む順がこうなる。
                        if (ICN2(1,II)==0.or.ICN2(5,II)==0) print*, 'ICN2_WARNING_pyr:', ICN2(:,II)
        
                    end if
    
                ELSE
                    read(n_unit,'()')  !テトラでもプリズムでもピラミッドでもないならスルー
    
                ENDIF
    
            END DO

            IIMX = II
        
            allocate(UVWK(3,KKMX))
    
            read(n_unit,*)AAmax
            read(n_unit,'()')
            DO II = 1,AAmax
                read(n_unit,'()')
            END DO
            DO KK = 1, KKMX
                read(n_unit,*)AA,UVWK(1,KK),UVWK(2,KK),UVWK(3,KK)
            END DO
                
        close(n_unit)

        if(.not.allocated(ICN)) allocate(ICN(6,IIMX))
        if(.not.allocated(CELL_TYPE)) allocate(CELL_TYPE(IIMX))
        ICN(:,:) = ICN2(:,:IIMX)
        CELL_TYPE(:) = CELL_TYPE2(:IIMX)

        if(.not.allocated(VELC)) allocate(VELC(3,IIMX))
            
        call point2cell(UVWK(:,:), VELC(:,:), ICN(:,:), CELL_TYPE(:))
            
    end subroutine read_INP

    subroutine read_FLD(FNAME)
        use mod_SctFldReader
        implicit none
        character(*), intent(in) :: FNAME
        integer unit, II,KK,KK_beg,KK_end,node_num
        double precision, allocatable :: UVWK(:,:)

  
        call open_readFLD(unit, FNAME)
            call read_Header_data(unit)
            call read_Main_data(unit)
        call close_fld(unit)
  
        if(.not.allocated(ICN)) then
            allocate(ICN(6,size(ietyp)), source=0)
            allocate(CELL_TYPE(size(ietyp)), source=0)
            allocate(CDN(3,NNODS))
  
            KK_beg = 1
  
            do II = 1, size(ietyp)
  
                node_num = ietyp(II)-30
                if (node_num==4) then
                    CELL_TYPE(II) = 0
                elseif(node_num==6) then
                    CELL_TYPE(II) = 1
                elseif(node_num==5) then
                    CELL_TYPE(II) = 2
                else
                    print*, 'ERROR_node_num', node_num
                end if
    
                KK_end = KK_beg + node_num - 1
    
                do KK = KK_beg, KK_end
                    ICN(KK-KK_beg+1, II) = ndno(KK) + 1
                end do
            
                KK_beg = KK_beg + node_num
        
            end do
        end if
  
        CDN(1,:) = CDN_X(:)
        CDN(2,:) = CDN_Y(:)
        CDN(3,:) = CDN_Z(:)

        allocate(UVWK(3, size(CDN_X)))
        UVWK(1,:) = VEL_X(:)
        UVWK(2,:) = VEL_Y(:)
        UVWK(3,:) = VEL_Z(:)

        if(.not.allocated(VELC)) allocate(VELC(3, size(ietyp)))

        call point2cell(UVWK(:,:), VELC(:,:), ICN(:,:), CELL_TYPE(:))
          
    end subroutine read_FLD
            
    !**************************************************************************************
    
    !***********************************************************************
    subroutine set_gravity_center !セル重心の算出
        integer II,IIMX
        integer ICN1,ICN2,ICN3,ICN4,ICN5,ICN6
            !=======================================================================
        IIMX = size(ICN(:,:), dim=2)    

        !$omp parallel do private(ICN1,ICN2,ICN3,ICN4,ICN5,ICN6)
        DO II = 1, IIMX
            IF (CELL_TYPE(II)==0) THEN
                ICN1 = ICN(1,II)
                ICN2 = ICN(2,II) 
                ICN3 = ICN(3,II)
                ICN4 = ICN(4,II)
            
                CENC(1,II) = 0.25d0*(CDN(1,ICN1)+CDN(1,ICN2)+CDN(1,ICN3)+CDN(1,ICN4))
                CENC(2,II) = 0.25d0*(CDN(2,ICN1)+CDN(2,ICN2)+CDN(2,ICN3)+CDN(2,ICN4))
                CENC(3,II) = 0.25d0*(CDN(3,ICN1)+CDN(3,ICN2)+CDN(3,ICN3)+CDN(3,ICN4))
            ELSE IF (CELL_TYPE(II)==1) THEN
                ICN1 = ICN(1,II)
                ICN2 = ICN(2,II) 
                ICN3 = ICN(3,II)
                ICN4 = ICN(4,II)
                ICN5 = ICN(5,II)
                ICN6 = ICN(6,II)
            
                CENC(1,II) = (CDN(1,ICN1)+CDN(1,ICN2)+CDN(1,ICN3)+CDN(1,ICN4)+CDN(1,ICN5)+CDN(1,ICN6))/6.0d0
                CENC(2,II) = (CDN(2,ICN1)+CDN(2,ICN2)+CDN(2,ICN3)+CDN(2,ICN4)+CDN(2,ICN5)+CDN(2,ICN6))/6.0d0
                CENC(3,II) = (CDN(3,ICN1)+CDN(3,ICN2)+CDN(3,ICN3)+CDN(3,ICN4)+CDN(3,ICN5)+CDN(3,ICN6))/6.0d0    
            ELSE IF (CELL_TYPE(II)==2) THEN
                ICN1 = ICN(1,II)
                ICN2 = ICN(2,II) 
                ICN3 = ICN(3,II)
                ICN4 = ICN(4,II)
                ICN5 = ICN(5,II)
            
                CENC(1,II) = 0.20d0*(CDN(1,ICN1)+CDN(1,ICN2)+CDN(1,ICN3)+CDN(1,ICN4)+CDN(1,ICN5))
                CENC(2,II) = 0.20d0*(CDN(2,ICN1)+CDN(2,ICN2)+CDN(2,ICN3)+CDN(2,ICN4)+CDN(2,ICN5))
                CENC(3,II) = 0.20d0*(CDN(3,ICN1)+CDN(3,ICN2)+CDN(3,ICN3)+CDN(3,ICN4)+CDN(3,ICN5))

            else
                print*, 'CELL_TYPE_ERROR', CELL_TYPE(II)
                stop
            END IF
            
            WIDC(II) = norm2(CDN(:,ICN2)-CDN(:,ICN1))
        
        END DO
        !$omp end parallel do 
    end subroutine set_gravity_center

    subroutine read_adjacency(path, success)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        logical, optional :: success
        integer II,NC,JB, n_unit, num_cells, NCMAX
        character(len_trim(path)+20) FNAME
                
        FNAME = trim(path)//adjacencyFileName
        if(present(success)) then
            inquire(file = FNAME, exist = success)
            if(.not.success) return
        end if

        print*, 'READ:', FNAME

        open(newunit=n_unit, FILE=FNAME, STATUS='OLD')
            read(n_unit,*) num_cells
            read(n_unit,*) NCMAX

            allocate(NEXT_CELL(NCMAX,num_cells),NUM_NC(num_cells))
            DO II = 1, num_cells
                read(n_unit,'(I5)',advance='no') NUM_NC(II)
                DO NC = 1, NUM_NC(II)
                    read(n_unit,'(I12)',advance='no') NEXT_CELL(NC,II)
                END DO
                read(n_unit,'()')
            END DO

            if(.not.allocated(NoB)) then

                allocate(NoB(num_cells), source=0)
                allocate(ICB(4,num_cells), source=0)

                DO II = 1, num_cells
                read(n_unit, fmt='(I4)', advance='no') NoB(II)  ! Number of Boundary
                do JB = 1, NoB(II)
                    read(n_unit, fmt='(I10)', advance='no') ICB(JB,II)
                end do
                read(n_unit,'()')  !改行
                END DO

            end if

        close(n_unit)

    end subroutine read_adjacency

    subroutine output_adjacency(path)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        integer II,NC,JB, n_unit, num_cells, NCMAX
        character(len_trim(path)+20) FNAME
                
        FNAME = trim(path)//adjacencyFileName

        print*, 'READ:', FNAME

        num_cells = size(NEXT_CELL(:,:), dim=2)
        NCMAX = size(NEXT_CELL(:,:), dim=1)

        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit,*) num_cells
            write(n_unit,*) NCMAX

            DO II = 1, num_cells
                write(n_unit,'(I5)',advance='no') NUM_NC(II)
                DO NC = 1, NUM_NC(II)
                    write(n_unit,'(I12)',advance='no') NEXT_CELL(NC,II)
                END DO
                write(n_unit,'()')
            END DO

            DO II = 1, num_cells
                write(n_unit, fmt='(I4)', advance='no') NoB(II)  ! Number of Boundary
                do JB = 1, NoB(II)
                    write(n_unit, fmt='(I10)', advance='no') ICB(JB,II)
                end do
                write(n_unit,'()')  !改行
            END DO

        close(n_unit)

    end subroutine output_adjacency

    subroutine read_boundaries(path)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(len_trim(path)+20) FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'READ:', FNAME
        open(newunit=n_unit, FILE=FNAME , STATUS='old')
            read(n_unit,*) JBMX
            allocate(NBN(3,JBMX), source=-1)
            do JB = 1, JBMX
                read(n_unit,*) NBN(1,JB), NBN(2,JB), NBN(3,JB)
            end do
        close(n_unit)
        
    end subroutine read_boundaries

    subroutine output_boundaries(path)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        integer JB, n_unit, JBMX
        character(len_trim(path)+20) FNAME

        FNAME = trim(path)//boundaryFileName
        print*, 'READ:', FNAME
        JBMX = size(NBN(:,:), dim=2)
        open(newunit=n_unit, FILE=FNAME , STATUS='replace')
            write(n_unit,*) JBMX
            do JB = 1, JBMX
                write(n_unit,*) NBN(1,JB), NBN(2,JB), NBN(3,JB)
            end do
        close(n_unit)
        
    end subroutine output_boundaries

    subroutine set_nodeID_onBoundFace(nodeID_onBoundFace)
        integer, intent(in) :: nodeID_onBoundFace(:,:)

        NBN = nodeID_onBoundFace
    end subroutine set_nodeID_onBoundFace

    subroutine set_adjacency(num_adjacent, adjacentCellID)
        integer, intent(in) :: num_adjacent(:), adjacentCellID(:,:)

        NUM_NC = num_adjacent
        NEXT_CELL = adjacentCellID
    end subroutine set_adjacency


    integer function nearest_cell(X) !最も近いセルNCNの探索
        double precision, intent(in) :: X(3)
        integer II, IIMX
        double precision, allocatable :: distance(:)

        IIMX = size(CENC, dim=2)
        allocate(distance(IIMX))

        !$omp parallel do
        DO II = 1,IIMX
                distance(II) = norm2(CENC(:,II) - X(:))
        END DO
        !$omp end parallel do 
        
        nearest_cell = minloc(distance, dim=1)   !最小値インデックス
        
    end function nearest_cell

    integer function nearer_cell(X, NCN)  !近セルの探索（隣接セルから）
        integer, intent(in) :: NCN
        double precision, intent(in) :: X(3)
        integer NC, IIaround, index_min
        double precision :: distancecheck(2)
        double precision, allocatable :: distance(:)

        nearer_cell = NCN
        allocate(distance(size(NEXT_CELL(:,:), dim=1)))
        distancecheck(1) = norm2(CENC(:,nearer_cell)-X(:))   !注目セル重心と粒子との距離
        
        check:DO
                distance(:) = 1.0d10     !初期値はなるべく大きくとる
        
                DO NC = 1, NUM_NC(nearer_cell)  !全隣接セルに対してループ。
                    IIaround = NEXT_CELL(NC, nearer_cell)       !現時点で近いとされるセルの隣接セルのひとつに注目
                    IF (IIaround > 0) then
                            distance(NC) = norm2(CENC(:,IIaround)-X(:))   !注目セル重心と粒子との距離を距離配列に代入
                    END IF
        
                END DO
        
                distancecheck(2) = minval(distance,dim=1)     !距離配列の最小値
        
                if(distancecheck(2) < distancecheck(1)) then !より近いセルの発見で条件満足
                    distancecheck(1) = distancecheck(2)    !最小値の更新
                    index_min = minloc(distance,dim=1)            !最小値のインデックス
                    nearer_cell = NEXT_CELL(index_min, nearer_cell)    !現時点で近いとされるセルの更新
                    if(nearer_cell==0) then
                            print*,'nearer_cell_error', nearer_cell, X(:)
                            return
                    end if

                else  !より近いセルを発見できなかった場合
        
                    exit check     !ループ脱出
        
                end if
        
        END DO check
        
    end function nearer_cell

    logical function nearcell_check(X, NCN)
        double precision, intent(in) :: X(3)
        integer, intent(in) :: NCN
        double precision :: distance

        distance = norm2(X(:)-CENC(:,NCN))

        if (distance < 1.0d1*WIDC(NCN)) then
            nearcell_check = .True.
        else
            nearcell_check = .False.
        end if


    end function nearcell_check
                     
    subroutine boundary_setting !全境界面に対して外向き法線ベクトルと重心を算出
        integer II, JJ, JB, IIMX, JBMX
        double precision :: a(3), b(3), r(3), norm, inner
        double precision, allocatable :: CENF_pre(:,:)

        if(.not.allocated(NoB)) return

        IIMX = size(ICN(:,:), dim=2)

        if(.not.allocated(CENF)) then
            JBMX = size(NBN(:,:), dim=2)
            allocate(CENF(3,JBMX), NVECF(3,JBMX))
        else
            allocate(CENF_pre, source=CENF)
        end if
        
        print*, 'SET:boundary'
        
        do II = 1, IIMX
            
            do JJ = 1, NoB(II)
                JB = ICB(JJ,II)
                
                CENF(:,JB) = (CDN(:,NBN(1,JB)) + CDN(:,NBN(2,JB)) + CDN(:,NBN(3,JB)))/3.0d0
            
                a(:) =  CDN(:,NBN(2,JB)) - CDN(:,NBN(1,JB))
                b(:) =  CDN(:,NBN(3,JB)) - CDN(:,NBN(1,JB))
            
                NVECF(1,JB) = a(2)*b(3) - a(3)*b(2)  !外積
                NVECF(2,JB) = a(3)*b(1) - a(1)*b(3)
                NVECF(3,JB) = a(1)*b(2) - a(2)*b(1)
            
                norm = norm2(NVECF(:,JB))
            
                r(:) = CENC(:,II) - CENF(:,JB)  !面重心からセル重心へのベクトル
            
                inner = sum(NVECF(:,JB)*r(:))
            
                if (inner > 0.0d0) norm = norm * (-1.0d0) !内積が正なら内向きなので、外に向けるべくノルムを負に
            
                NVECF(:,JB) = NVECF(:,JB) / norm !ノルムで割り算して単位ベクトルに
        
            end do
        
        end do

        if(allocated(CENF_pre)) MOVF = CENF - CENF_pre  !自動割付

    end subroutine boundary_setting

    subroutine point2cell(pointdata, celldata, cell2node, celltype)
        implicit none
        double precision, intent(in) :: pointdata(:,:)
        double precision, intent(inout) :: celldata(:,:)
        integer, intent(in) :: cell2node(:,:), celltype(:)
        integer II
  
        do II = 1, size(celldata, dim=2)   !点データをセルデータに変換
              
            IF (celltype(II)==0) THEN
  
                celldata(:,II) = 0.25d0*(pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) &
                    + pointdata(:, cell2node(3,II)) + pointdata(:, cell2node(4,II)))
  
            ELSE IF (celltype(II)==1) THEN
  
                celldata(:,II) = (pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) + pointdata(:, cell2node(3,II)) &
                    + pointdata(:, cell2node(4,II)) + pointdata(:, cell2node(5,II)) + pointdata(:, cell2node(6,II))) / 6.0d0
  
            ELSE IF (celltype(II)==2) THEN
  
                celldata(:,II) = 0.20d0*(pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) &
                    + pointdata(:, cell2node(3,II)) + pointdata(:, cell2node(4,II)) + pointdata(:, cell2node(5,II)))
  
            END IF
  
        end do
  
    end subroutine point2cell

    integer function get_mesh_info(name)
        character(*), intent(in) :: name
        
        select case(name)
            case('node')
                get_mesh_info = size(CDN(:,:), dim = 2)

            case('cell')
                get_mesh_info = size(ICN(:,:), dim = 2)

            case('tetra')
                get_mesh_info = count(CELL_TYPE == 0)

            case('prism')
                get_mesh_info = count(CELL_TYPE == 1)

            case('pyramid')
                get_mesh_info = count(CELL_TYPE == 2)

            case default
                get_mesh_info = 0

        end select

    end function get_mesh_info
      
    subroutine deallocation_flow

        if(allocated(CDN)) deallocate(CDN, ICN, CELL_TYPE, VELC)
        if(allocated(NEXT_CELL)) deallocate(NEXT_CELL, NUM_NC)
        if(allocated(NoB)) deallocate(NoB, ICB)
        if(allocated(NBN)) deallocate(NBN)
        if(allocated(CENF)) deallocate(CENF, NVECF)
        if(allocated(CENC)) deallocate(CENC, WIDC)

    end subroutine deallocation_flow
    
end module unstructuredGrid_mod