module flow_field
    implicit none
    integer, private :: NCMAX                                   !隣接要素数の最大値
    integer INTERVAL_FLOW                           !気流データ出力間隔

    character PATH_AIR*99, HEAD_AIR*20, FNAME_FMT*30 !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    character(3) FILE_TYPE  !ファイル形式（VTK,INP）
    integer, private :: FNAME_DIGITS !ファイル名の整数部桁数

    integer, allocatable :: ICN(:,:)                !要素所有節点ID
    integer, allocatable :: NoB(:)                  !要素所有境界面の数
    integer, allocatable :: ICB(:,:)                !要素所有境界面ID
    integer, allocatable :: NBN(:,:)                !境界面所有節点ID
    
    integer, allocatable :: CELL_TYPE(:)   !要素タイプ（テトラ0、プリズム1、ピラミッド2）
    integer, allocatable, private :: NUM_NC(:)               !隣接要素数
    integer, allocatable, private :: NEXT_CELL(:,:)          !隣接要素ID

    double precision, allocatable :: CDN(:,:)       !節点座標
    double precision MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限
    double precision, allocatable :: VELC(:,:)      !要素流速
    double precision, allocatable, private :: CENC(:,:)      !要素重心
    double precision, allocatable, private :: WIDC(:)        !要素の1辺長さ
    double precision, allocatable :: CENF(:,:,:)    !面重心
    double precision, allocatable :: NVECF(:,:)     !面法線ベクトル

    contains
    !***********************************************************************
    subroutine set_FILE_TYPE
        integer i

        i = index(FNAME_FMT, '0')
        HEAD_AIR = FNAME_FMT(: i-1)             !ファイル名の接頭文字(最初のゼロの手前まで)
        FNAME_DIGITS = index(FNAME_FMT, '.') - i   !ファイル名の整数部桁数(最初のゼロの位置からドットまでの文字数)

        if(index(FNAME_FMT, '.vtk') > 0) then
            FILE_TYPE = 'VTK'

        else if(index(FNAME_FMT, '.inp') > 0) then
            FILE_TYPE = 'INP'

        else if(index(FNAME_FMT, '.fld') > 0) then
            FILE_TYPE = 'FLD'

        else
            print*, 'FNAME_FMT NG:', FNAME_FMT
            stop
        end if

        print*, 'FILE_TYPE: ', FILE_TYPE, ' ', trim(FNAME_FMT)

    end subroutine set_FILE_TYPE

    character(4) function get_digits_format()

        write(get_digits_format,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS

    end function get_digits_format

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

    subroutine read_nextcell
        implicit none
        integer II,NC,JB, n_unit, num_cells, JBMX
        character(50) FNAME
                
        FNAME = trim(PATH_AIR)//'nextcell.txt'
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

            allocate(NoB(num_cells), source=0)
            allocate(ICB(4,num_cells), source=0)
            allocate(CENC(3,num_cells), WIDC(num_cells))

            DO II = 1, num_cells
            read(n_unit, fmt='(I4)', advance='no') NoB(II)  ! Number of Boundary
            do JB = 1, NoB(II)
                read(n_unit, fmt='(I10)', advance='no') ICB(JB,II)
            end do
            read(n_unit,'()')  !改行
            END DO

        close(n_unit)


        FNAME = trim(PATH_AIR)//'boundaries.txt'
        print*, 'READ:', FNAME
        open(newunit=n_unit, FILE=FNAME , STATUS='old')
            read(n_unit,*) JBMX
            allocate(NBN(3,JBMX), source=-1)
            do JB = 1, JBMX
                read(n_unit,*) NBN(1,JB), NBN(2,JB), NBN(3,JB)
            end do
        close(n_unit)
        allocate(CENF(3,JBMX,2), NVECF(3,JBMX))
        
    end subroutine read_nextcell

    integer function nearest_cell(X) !最も近いセルNCNの探索
        double precision, intent(in) :: X(3)
        integer II, IIMX
        double precision, allocatable :: distance(:)

        IIMX = size(CENC, dim=2)

        allocate(distance(IIMX))
        !↓↓↓↓　一番近いセル中心の探索
        !$omp parallel do
        DO II = 1,IIMX
                distance(II) = norm2(CENC(:,II) - X(:))
        END DO
        !$omp end parallel do 
        !↑↑↑↑
        
        nearest_cell = minloc(distance, dim=1)   !最小値インデックス
        
    end function nearest_cell

    integer function nearer_cell(X, NCN)  !近セルの探索（隣接セルから）
        integer, intent(in) :: NCN
        double precision, intent(in) :: X(3)
        integer NC, IIaround, index_min
        double precision :: distancecheck(2)
        double precision, allocatable :: distance(:)

        nearer_cell = NCN
        allocate(distance(NCMAX))
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
                     
    subroutine boundary_set !全境界面に対して法線ベクトルと重心を算出
        integer II, JJ, JB, IIMX
        double precision :: a(3), b(3), r(3), norm, inner

        IIMX = size(ICN(:,:), dim=2)
        
        print*, 'SET:boundary'
        
        do II = 1, IIMX
            
            do JJ = 1, NoB(II)
                JB = ICB(JJ,II)
                
                CENF(:,JB,2) = (CDN(:,NBN(1,JB)) + CDN(:,NBN(2,JB)) + CDN(:,NBN(3,JB)))/3.0d0
            
                a(:) =  CDN(:,NBN(2,JB)) - CDN(:,NBN(1,JB))
                b(:) =  CDN(:,NBN(3,JB)) - CDN(:,NBN(1,JB))
            
                NVECF(1,JB) = a(2)*b(3) - a(3)*b(2)  !外積
                NVECF(2,JB) = a(3)*b(1) - a(1)*b(3)
                NVECF(3,JB) = a(1)*b(2) - a(2)*b(1)
            
                norm = norm2(NVECF(:,JB))
            
                r(:) = CENC(:,II) - CENF(:,JB,2)  !面重心からセル重心へのベクトル
            
                inner = sum(NVECF(:,JB)*r(:))
            
                if (inner > 0.0d0) norm = norm * (-1.0d0) !内積が正なら内向きなので、外に向けるべくノルムを負に
            
                NVECF(:,JB) = NVECF(:,JB) / norm !ノルムで割り算して単位ベクトルに
        
            end do
        
        end do
        
        !----------------------------------------------------------------------------------
    end subroutine boundary_set
        !----------------------------------------------------------------------------------

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

        deallocate(CDN, ICN, CELL_TYPE, VELC)
        if(allocated(NEXT_CELL)) deallocate(NEXT_CELL, NUM_NC)
        deallocate(NoB, ICB)
        if(allocated(NBN)) deallocate(NBN)
        if(allocated(CENF)) deallocate(CENF, NVECF)
        if(allocated(CENC)) deallocate(CENC, WIDC)

    end subroutine deallocation_flow
    
end module flow_field