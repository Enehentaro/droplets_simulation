module flow_field
    implicit none
    integer, private :: IIMX, KKMX, JBMX            !全要素数、全節点数、全境界面数
    integer NCMAX                                   !隣接要素数の最大値
    integer INTERVAL_FLOW                           !気流データ出力間隔

    character PATH_AIR*99, HEAD_AIR*10, FNAME_FMT*20 !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    character(3) FILE_TYPE  !ファイル形式（VTK,INP）
    integer FNAME_DIGITS !ファイル名の整数部桁数

    integer, allocatable :: ICN(:,:)                !要素所有節点ID
    integer, allocatable :: NoB(:)                  !要素所有境界面の数
    integer, allocatable :: ICB(:,:)                !要素所有境界面ID
    integer, allocatable :: NBN(:,:)                !境界面所有節点ID
    
    integer, allocatable :: CELL_TYPE(:)   !要素タイプ（テトラ0、プリズム1、ピラミッド2）
    integer, allocatable :: NUM_NC(:)               !隣接要素数
    integer, allocatable :: NEXT_CELL(:,:)          !隣接要素ID

    double precision, allocatable :: CDN(:,:)       !節点座標
    double precision MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限
    double precision, allocatable :: VELC(:,:)      !要素流速
    double precision, allocatable :: CENC(:,:)      !要素重心
    double precision, allocatable :: WIDC(:)        !要素の1辺長さ
    double precision, allocatable :: CENF(:,:,:)    !面重心
    double precision, allocatable :: NVECF(:,:)     !面法線ベクトル

    contains
    !***********************************************************************
    subroutine set_FILE_TYPE
        integer i

        i = index(FNAME_FMT, '0')
        HEAD_AIR = FNAME_FMT(: i-1)             !ファイル名の接頭文字
        FNAME_DIGITS = index(FNAME_FMT, '.') - i   !ファイル名の整数部桁数

        if(index(FNAME_FMT, '.vtk') > 0) then
            FILE_TYPE = 'VTK'

        else if(index(FNAME_FMT, '.inp') > 0) then
            FILE_TYPE = 'INP'

        else
            print*, 'FNAME_FMT NG:', FNAME_FMT
            stop
        end if

    end subroutine set_FILE_TYPE


    subroutine read_VTK(FNAME, pointdata)
        character(*), intent(in) :: FNAME
        logical, optional :: pointdata
        integer II,KK,IIH, n_unit
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
        integer :: IITETMX, IIPRSMX, IIPYRMX
        character(6) cellshape
        double precision, allocatable :: UVWK(:,:)

        print*, 'READ_INP:', trim(FNAME)
            
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
            read(n_unit,*)KKMX,IIMX2
            print*,'KKMX,IIMX2=',KKMX,IIMX2
            
            if(.not.allocated(CDN)) allocate(CDN(3,KKMX), source=0.0d0)
            if(.not.allocated(ICN)) allocate(ICN(6,IIMX2), source=0)
            if(.not.allocated(CELL_TYPE)) allocate(CELL_TYPE(IIMX2), source=0)
                
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
                        CELL_TYPE(II) = 0
                        IITETMX = IITETMX +1
                        read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II)    
                        if (ICN(1,II)==0.or.ICN(4,II)==0) print*, 'ICN_WARNING_tet:', ICN(:,II)
                    
                    ELSE IF(cellshape=='prism') THEN
                        CELL_TYPE(II) = 1
                        IIPRSMX = IIPRSMX +1
                        read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II),ICN(5,II),ICN(6,II)
                        if (ICN(1,II)==0.or.ICN(6,II)==0) print*, 'ICN_WARNING_prism:', ICN(:,II)
        
                    ELSE IF(cellshape=='pyr') THEN
                        CELL_TYPE(II) = 2
                        IIPYRMX = IIPYRMX +1
                        read(n_unit,*)ICN(5,II),ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II) !INPは最初が山頂点であり、VTKでは最後が山頂点のため、読み込む順がこうなる。
                        if (ICN(1,II)==0.or.ICN(5,II)==0) print*, 'ICN_WARNING_pyr:', ICN(:,II)
        
                    end if
    
                ELSE
                    read(n_unit,'()')  !テトラでもプリズムでもピラミッドでもないならスルー
    
                ENDIF
    
            END DO
    
            if (II == IIMX) then
                print*, 'IIMX=', IIMX
                print*, 'Tetra,Prism,Pyramid=', IITETMX, IIPRSMX, IIPYRMX
            else
                print*, 'IIMX_mismatch', II, '/', IIMX
                stop
            end if
        
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
            
        call point2cell(UVWK(:,:), VELC(:,:), ICN(:,:), CELL_TYPE(:))
            
    end subroutine read_INP
            
    !**************************************************************************************
    
    !***********************************************************************
    subroutine set_gravity_center !セル重心の算出
        integer II
        integer ICN1,ICN2,ICN3,ICN4,ICN5,ICN6
            !=======================================================================
            
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


            !***********************************************************************
    subroutine read_nextcell
        integer II,NC,JB, n_unit
        character(50) FNAME
        !=======================================================================
                
        FNAME = trim(PATH_AIR)//'nextcell.txt'
        print*, 'READ:', FNAME

        open(newunit=n_unit,FILE=FNAME , STATUS='OLD')
            read(n_unit,'(2(I12,2X))') IIMX

            read(n_unit,*) NCMAX

            allocate(NEXT_CELL(NCMAX,IIMX),NUM_NC(IIMX))
            DO II = 1, IIMX
            read(n_unit,'(I5)',advance='no') NUM_NC(II)
            DO NC = 1, NUM_NC(II)
                read(n_unit,'(I12)',advance='no') NEXT_CELL(NC,II)
            END DO
            read(n_unit,'()')
            END DO

            allocate(NoB(IIMX), source=0)
            allocate(ICB(4,IIMX), source=0)
            allocate(CENC(3,IIMX), WIDC(IIMX))

            ! DO II = 1, IIMX
            !   read(n_unit,*) NoB(II)  ! Number of Boundary
            ! END DO
            DO II = 1, IIMX
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
        !*************************************************************************************
        !----------------------------------------------------------------------------------                        
    subroutine boundary_set !全境界面に対して法線ベクトルと重心を算出
        integer II, JJ, JB
        double precision :: a(3), b(3), r(3), norm, inner
        
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

    integer function get_num_nodes()

        get_num_nodes = KKMX

    end function get_num_nodes

      
    subroutine deallocation_flow

        deallocate(CDN, ICN, CELL_TYPE, VELC)
        if(allocated(NEXT_CELL)) deallocate(NEXT_CELL, NUM_NC)
        deallocate(NoB, ICB)
        if(allocated(NBN)) deallocate(NBN)
        if(allocated(CENF)) deallocate(CENF, NVECF)
        if(allocated(CENC)) deallocate(CENC, WIDC)

    end subroutine deallocation_flow
    
end module flow_field