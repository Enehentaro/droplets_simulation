module unstructuredGrid_mod
    implicit none
    character(3), private :: FILE_TYPE  !ファイル形式

    type node_t
        real coordinate(3)
    end type node_t

    type boundFace_t
        integer nodeID(3)
        real center(3), normalVector(3), moveVector(3)
    end type boundFace_t

    type cell_t
        character(5) typeName
        integer, allocatable :: nodeID(:), boundFaceID(:), adjacentCellID(:)
        real center(3), flowVelocity(3), width
    end type cell_t

    type(node_t), allocatable :: NODEs(:)
    type(boundFace_t), allocatable :: BoundFACEs(:)
    type(cell_t), allocatable :: CELLs(:)

    ! integer, allocatable :: ICN(:,:)                !要素所有節点ID
    ! integer, allocatable :: NoB(:)                  !要素所有境界面の数
    ! integer, allocatable :: ICB(:,:)                !要素所有境界面ID
    ! integer, allocatable :: NBN(:,:)                !境界面所有節点ID
    
    ! integer, allocatable :: CELL_TYPE(:)   !要素タイプ（テトラ0、プリズム1、ピラミッド2）
    ! integer, allocatable, private :: NUM_NC(:)               !隣接要素数
    ! integer, allocatable, private :: NEXT_CELL(:,:)          !隣接要素ID

    ! double precision, allocatable :: CDN(:,:)       !節点座標
    ! double precision, allocatable :: VELC(:,:)      !要素流速
    ! double precision, allocatable :: CENC(:,:)      !要素重心
    ! double precision, allocatable, private :: WIDC(:)        !要素の1辺長さ
    ! double precision, allocatable :: CENF(:,:)    !面重心
    ! double precision, allocatable :: MOVF(:,:)    !面重心移動量
    ! double precision, allocatable :: NVECF(:,:)     !面法線ベクトル

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
            
    end subroutine read_unstructuredGrid_byNAME

    subroutine read_unstructuredGrid_byNumber(path_and_head, digits_fmt, FNUM)
        character(*), intent(in) :: path_and_head
        integer, intent(in) :: FNUM
        character(99) :: FNAME
        character(4) digits_fmt

        select case(FILE_TYPE)
            case('VTK')

                write(FNAME,'("'//trim(path_and_head)//'",'//digits_fmt//',".vtk")') FNUM
                call read_VTK(FNAME)

            case('INP')

                if(FNUM==0) then
                    write(FNAME,'("'//trim(path_and_head)//'",'//digits_fmt//',".inp")') 1
                else
                    write(FNAME,'("'//trim(path_and_head)//'",'//digits_fmt//',".inp")') FNUM
                end if

                call read_INP(FNAME)   !INPを読み込む(SHARP用)

            case('FLD')
                if(.not.allocated(CELLs).and.FNUM > 0) then
                    write(FNAME,'("'//trim(path_and_head)//'", i0, ".fld")') 0
                    call read_FLD(FNAME)
                end if

                write(FNAME,'("'//trim(path_and_head)//'", i0, ".fld")') FNUM
                call read_FLD(FNAME)

            case default
                print*,'FILE_TYPE NG:', FILE_TYPE
                STOP
                
        end select
            
    end subroutine read_unstructuredGrid_byNumber

    subroutine read_VTK(FNAME, pointdata)
        use vtkMesh_operator_m
        character(*), intent(in) :: FNAME
        logical, optional :: pointdata
        integer II,KK,IIH, KKMX, IIMX
        real, allocatable :: UVWK(:,:)

        call read_VTK_mesh(FNAME)

        KKMX = size(node_array)
        if(.not.allocated(NODEs)) allocate(NODEs(KKMX))
        do KK = 1, KKMX
            NODEs(KK)%coordinate(:) = node_array(KK-1)%coordinate(:)
        end do
        
        IIMX = size(cell_array)
        if(.not.allocated(CELLs)) allocate(CELLs(IIMX))
        do II = 1, IIMX
            IIH = size(cell_array(II-1)%nodeID)
            CELLs(II)%nodeID = cell_array(II-1)%nodeID(1:IIH) + 1
            select case(cell_array(II-1)%n_TYPE)
                case(10)
                    CELLs(II)%typeName = 'tetra'
                case(13)
                    CELLs(II)%typeName = 'prism'
                case(14)
                    CELLs(II)%typeName = 'pyrmd'
            end select
            CELLs(II)%flowVelocity(:) = cell_array(II-1)%vector(:)
        end do

        ! print*, NODEs(KKMX)%coordinate(:)
        ! print*, CELLs(IIMX)%flowVelocity(:)

        if(present(pointdata)) then
            if(pointdata) then
                allocate(UVWK(3,KKMX))

                call point2cellVelocity(UVWK)

                return
            end if
        end if
            
    end subroutine read_VTK

    subroutine read_INP(FNAME)
        !  INPファイルを読み込み、節点データを要素データに変換する
        character(*), intent(in) :: FNAME
        INTEGER II,II2,KK,AAmax, IIMX2, AA, n_unit
        integer KKMX, IIMX
        character(6) cellshape
        real, allocatable :: UVWK(:,:)
        integer, allocatable :: ICN2(:,:)
        integer, allocatable :: CELL_TYPE2(:)

        print*, 'READ_INP:', trim(FNAME)
            
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
            read(n_unit,*)KKMX,IIMX2
            print*,'KKMX,IIMX2=',KKMX,IIMX2
            
            if(.not.allocated(NODEs)) allocate(NODEs(KKMX))
            allocate(ICN2(6,IIMX2), source=0)
            allocate(CELL_TYPE2(IIMX2), source=0)
                
            DO KK = 1, KKMX
                read(n_unit,*)AA, NODEs(KK)%coordinate(:)
            END DO
                
            II = 0
            DO II2 = 1, IIMX2
    
                read(n_unit,fmt='(I10)',advance='no')AA  !ここはセル番号なので無視
                read(n_unit,fmt='(I6)',advance='no')AA  !ここはなんかしらん（だいたいゼロ）
                read(n_unit,fmt='(A6)',advance='no')cellshape
    
                cellshape = adjustl(cellshape)    !左詰め
                IF ((cellshape=='tet').or.(cellshape=='prism').or.(cellshape=='pyr')) THEN
                    II = II +1
        
                    if (cellshape=='tet') then
                        CELL_TYPE2(II) = 0
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II)    
                        if (ICN2(1,II)==0.or.ICN2(4,II)==0) print*, 'ICN2_WARNING_tet:', ICN2(:,II)
                    
                    ELSE IF(cellshape=='prism') THEN
                        CELL_TYPE2(II) = 1
                        read(n_unit,*)ICN2(1,II),ICN2(2,II),ICN2(3,II),ICN2(4,II),ICN2(5,II),ICN2(6,II)
                        if (ICN2(1,II)==0.or.ICN2(6,II)==0) print*, 'ICN2_WARNING_prism:', ICN2(:,II)
        
                    ELSE IF(cellshape=='pyr') THEN
                        CELL_TYPE2(II) = 2
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
                read(n_unit,*)AA, UVWK(:,KK)
            END DO
                
        close(n_unit)

        if(.not.allocated(CELLs)) allocate(CELLs(IIMX))
        do II = 1, IIMX
            select case(CELL_TYPE2(II))
                case(0)
                    CELLs(II)%nodeID = ICN2(1:4, II)
                    CELLs(II)%typeName = 'tetra'
                case(1)
                    CELLs(II)%nodeID = ICN2(1:6, II)
                    CELLs(II)%typeName = 'prism'
                case(2)
                    CELLs(II)%nodeID = ICN2(1:5, II)
                    CELLs(II)%typeName = 'pyrmd'
            end select
        end do
            
        call point2cellVelocity(UVWK)
            
    end subroutine read_INP

    subroutine read_FLD(FNAME)
        use mod_SctFldReader
        implicit none
        character(*), intent(in) :: FNAME
        integer unit, II,KK,KK_beg,KK_end, num_node, IIMX
        real, allocatable :: UVWK(:,:)

        print*, 'readFLD : ', trim(FNAME)

  
        call open_readFLD(unit, FNAME)
            call read_Header_data(unit)
            call read_Main_data(unit)
        call close_fld(unit)
  
        if(.not.allocated(CELLs)) then
            allocate(CELLs(size(ietyp)))
            allocate(NODEs(NNODS))
  
            KK_beg = 1
            IIMX = size(ietyp)
  
            do II = 1, IIMX
                if(mod(II, 10000) == 0) print*, II, '/', IIMX
  
                num_node = ietyp(II)-30
                select case(num_node)
                    case(4)
                        CELLs(II)%typeName = 'tetra'
                    case(6)
                        CELLs(II)%typeName = 'prism'
                    case(5)
                        CELLs(II)%typeName = 'pyrmd'
                end select

                allocate(CELLs(II)%nodeID(num_node))
    
                KK_end = KK_beg + num_node - 1
    
                do KK = KK_beg, KK_end
                    CELLs(II)%nodeID(KK-KK_beg+1) = ndno(KK) + 1
                    ! ICN(KK-KK_beg+1, II) = ndno(KK) + 1
                end do
            
                KK_beg = KK_beg + num_node
        
            end do
        end if
  
        NODEs(:)%coordinate(1) = real(CDN_X(:))
        NODEs(:)%coordinate(2) = real(CDN_Y(:))
        NODEs(:)%coordinate(3) = real(CDN_Z(:))

        allocate(UVWK(3, size(CDN_X)))
        UVWK(1,:) = real(VEL_X(:))
        UVWK(2,:) = real(VEL_Y(:))
        UVWK(3,:) = real(VEL_Z(:))

        call point2cellVelocity(UVWK)
          
    end subroutine read_FLD
            
    !**************************************************************************************
    
    !***********************************************************************
    subroutine set_gravity_center !セル重心の算出
        integer II,IIMX, n, num_node, ID

        IIMX = size(CELLs)

        !$omp parallel do private(num_node, ID)
        DO II = 1, IIMX
            CELLs(II)%center(:) = 0.0
            num_node = size(CELLs(II)%nodeID)
            do n = 1, num_node
                ID = CELLs(II)%nodeID(n)
                CELLs(II)%center(:) = CELLs(II)%center(:) + NODEs(ID)%coordinate(:)
            end do
            CELLs(II)%center(:) = CELLs(II)%center(:) / num_node
            
            CELLs(II)%width = norm2( &
                NODEs(CELLs(II)%nodeID(2))%coordinate(:) - NODEs(CELLs(II)%nodeID(1))%coordinate(:))
        
        END DO
        !$omp end parallel do 
    end subroutine set_gravity_center

    subroutine read_adjacency(path, success)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        logical, optional :: success
        integer II,NA,JB, n_unit, num_cells, num_adj, num_BF, NCMAX
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

            ! allocate(NEXT_CELL(NCMAX,num_cells),NUM_NC(num_cells))
            DO II = 1, num_cells
                read(n_unit,'(I5)',advance='no') num_adj
                allocate(CELLs(II)%adjacentCellID(num_adj))
                DO NA = 1, num_adj
                    read(n_unit,'(I12)',advance='no') CELLs(II)%adjacentCellID(NA)
                END DO
                read(n_unit,'()')
            END DO

            DO II = 1, num_cells
                read(n_unit, fmt='(I4)', advance='no') num_BF  ! Number of Boundary
                allocate(CELLs(II)%boundFaceID(num_BF))
                do JB = 1, num_BF
                    read(n_unit, fmt='(I10)', advance='no') CELLs(II)%boundFaceID(JB)
                end do
                read(n_unit,'()')  !改行
            END DO


        close(n_unit)

    end subroutine read_adjacency

    subroutine output_adjacency(path)
        use filename_mod
        implicit none
        character(*), intent(in) :: path
        integer II,NA,JB, n_unit, num_cells, NCMAX
        character(len_trim(path)+20) FNAME
                
        FNAME = trim(path)//adjacencyFileName

        print*, 'OUTPUT:', FNAME

        num_cells = size(CELLs)
        NCMAX = 5   !size(NEXT_CELL(:,:), dim=1)

        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit,*) num_cells
            write(n_unit,*) NCMAX

            DO II = 1, num_cells
                write(n_unit,'(I5)',advance='no') size(CELLs(II)%adjacentCellID)
                DO NA = 1, size(CELLs(II)%adjacentCellID)
                    write(n_unit,'(I12)',advance='no') CELLs(II)%adjacentCellID(NA)
                END DO
                write(n_unit,'()')
            END DO

            DO II = 1, num_cells
                write(n_unit, fmt='(I4)', advance='no') size(CELLs(II)%boundFaceID)  ! Number of Boundary
                do JB = 1, size(CELLs(II)%boundFaceID)
                    write(n_unit, fmt='(I10)', advance='no') CELLs(II)%boundFaceID(JB)
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
            allocate(BoundFACEs(JBMX))
            do JB = 1, JBMX
                read(n_unit,*) BoundFACEs(JB)%nodeID(:)
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
        print*, 'OUTPUT:', FNAME
        JBMX = size(BoundFACEs)
        open(newunit=n_unit, FILE=FNAME , STATUS='replace')
            write(n_unit,*) JBMX
            do JB = 1, JBMX
                write(n_unit,*) BoundFACEs(JB)%nodeID(:)
            end do
        close(n_unit)
        
    end subroutine output_boundaries

    ! subroutine set_nodeID_onBoundFace(nodeID_onBoundFace)
    !     integer, intent(in) :: nodeID_onBoundFace(:,:)
    !     integer i, num_BF

    !     num_BF = size(nodeID_onBoundFace, dim=2)
    !     allocate(BoundFACEs(num_BF))

    !     do i = 1, num_BF
    !         BoundFACEs(i)%nodeID = nodeID_onBoundFace(:,i)
    !     end do
    ! end subroutine set_nodeID_onBoundFace

    ! subroutine set_adjacency(num_adjacent, adjacentCellID)
    !     integer, intent(in) :: num_adjacent(:), adjacentCellID(:,:)
    !     integer i, num_cell

    !     num_cell = size(CELLs)
    !     do i = 1, num_cell
    !         CELLs(i)%adjacentCellID = adjacentCellID(1:num_adjacent(i), i)
    !     end do

    ! end subroutine set_adjacency

    integer function nearest_cell(X) !最も近いセルNCNの探索
        real, intent(in) :: X(3)
        integer II, IIMX
        real, allocatable :: distance(:)

        IIMX = size(CELLs)
        allocate(distance(IIMX))

        !$omp parallel do
        DO II = 1,IIMX
            distance(II) = norm2(CELLs(II)%center(:) - X(:))
        END DO
        !$omp end parallel do 
        
        nearest_cell = minloc(distance, dim=1)   !最小値インデックス
        
    end function nearest_cell

    integer function nearer_cell(X, NCN)  !近セルの探索（隣接セルから）
        integer, intent(in) :: NCN
        real, intent(in) :: X(3)
        integer NA, featuredCELL, adjacentCELL
        real distance, distance_min
        logical update

        nearer_cell = NCN
        distance_min = norm2(CELLs(nearer_cell)%center(:) - X(:))   !注目セル重心と粒子との距離
        update = .true.
        do while(update)    !更新が起こり続ける限り繰り返し
            update = .false.
            featuredCELL = nearer_cell

            checkAdjacent : do NA = 1, size(CELLs(featuredCELL)%adjacentCellID)  !全隣接セルに対してループ。

                adjacentCELL = CELLs(featuredCELL)%adjacentCellID(NA)  !注目セルの隣接セルのひとつに注目
                if (adjacentCELL <= 0) cycle checkAdjacent

                distance = norm2(CELLs(adjacentCELL)%center(:) - X(:))   !注目セル重心と粒子との距離
                if(distance < distance_min) then
                    nearer_cell = adjacentCELL
                    distance_min = distance
                    update = .true.

                end if

            end do checkAdjacent

        end do
        
        ! check:DO
        !     distance(:) = 1.0d10     !初期値はなるべく大きくとる
    
        !     DO NA = 1, size(CELLs(nearer_cell)%adjacentCellID)  !全隣接セルに対してループ。
        !         IIaround = CELLs(nearer_cell)%adjacentCellID(NA)  !現時点で近いとされるセルの隣接セルのひとつに注目
        !         IF (IIaround > 0) distance(NA) = norm2(CENC(:,IIaround)-X(:))   !注目セル重心と粒子との距離を距離配列に代入
        !     END DO
    
        !     distancecheck(2) = minval(distance,dim=1)     !距離配列の最小値
    
        !     if(distancecheck(2) < distancecheck(1)) then !より近いセルの発見で条件満足
        !         distancecheck(1) = distancecheck(2)    !最小値の更新
        !         index_min = minloc(distance,dim=1)            !最小値のインデックス
        !         nearer_cell = NEXT_CELL(index_min, nearer_cell)    !現時点で近いとされるセルの更新
        !         if(nearer_cell==0) then
        !             print*,'nearer_cell_error', nearer_cell, X(:)
        !             return
        !         end if

        !     else  !より近いセルを発見できなかった場合
    
        !         exit check     !ループ脱出
    
        !     end if
        
        ! END DO check
        
    end function nearer_cell

    logical function nearcell_check(X, NCN)
        real, intent(in) :: X(3)
        integer, intent(in) :: NCN
        real :: distance

        distance = norm2(X(:) - CELLs(NCN)%center(:))

        if (distance < 1.0d1*CELLs(NCN)%width) then
            nearcell_check = .True.
        else
            nearcell_check = .False.
            print*, 'false', distance, CELLs(NCN)%width
        end if

    end function nearcell_check
                     
    subroutine boundary_setting(first) !全境界面に対して外向き法線ベクトルと重心を算出
        use vector_m
        logical, intent(in) :: first
        integer II, JJ, JB, IIMX, JBMX, nodeID(3)
        real :: a(3), b(3), r(3), norm, inner
        type(boundFace_t), allocatable :: BoundFACEs_pre(:)

        if(.not.allocated(BoundFACEs)) return

        IIMX = size(CELLs)

        if(.not.first) BoundFACEs_pre = BoundFACEs
        
        print*, 'SET:boundary'
        
        do II = 1, IIMX
            
            do JJ = 1, size(CELLs(II)%boundFaceID)
                JB = CELLs(II)%boundFaceID(JJ)
                nodeID(:) = BoundFACEs(JB)%nodeID(:)
                
                BoundFACEs(JB)%center(:) = ( NODEs(nodeID(1))%coordinate(:) &
                                            + NODEs(nodeID(2))%coordinate(:) &
                                            + NODEs(nodeID(3))%coordinate(:) ) / 3.0
            
                a(:) =  NODEs(nodeID(2))%coordinate(:) - NODEs(nodeID(1))%coordinate(:)
                b(:) =  NODEs(nodeID(3))%coordinate(:) - NODEs(nodeID(1))%coordinate(:)
            
                BoundFACEs(JB)%normalVector(1) = a(2)*b(3) - a(3)*b(2)  !外積
                BoundFACEs(JB)%normalVector(2) = a(3)*b(1) - a(1)*b(3)
                BoundFACEs(JB)%normalVector(3) = a(1)*b(2) - a(2)*b(1)
            
                norm = norm2(BoundFACEs(JB)%normalVector(:))
            
                r(:) = CELLs(II)%center(:) - BoundFACEs(JB)%center(:)  !面重心からセル重心へのベクトル
            
                inner = dot_product(BoundFACEs(JB)%normalVector(:), r(:))
            
                if(inner > 0.0) norm = norm * (-1.0) !内積が正なら内向きなので、外に向けるべくノルムを負に
            
                BoundFACEs(JB)%normalVector(:) = BoundFACEs(JB)%normalVector(:) / norm !ノルムで割り算して単位ベクトルに
        
                ! print*,'center:',BoundFACEs(JB)%center(:)
                ! print*,'n_vector:',BoundFACEs(JB)%normalVector(:)
            end do
        
        end do

        if(.not.first) then
            JBMX = size(BoundFACEs)
            do JB = 1, JBMX
                BoundFACEs(JB)%moveVector(:) = BoundFACEs(JB)%center(:) - BoundFACEs_pre(JB)%center(:)
            end do
        end if

    end subroutine boundary_setting

    subroutine point2cellVelocity(pointVector)
        real, intent(in) :: pointVector(:,:)
        integer II, IIMX, n, ID, num_node

        IIMX = size(CELLs)
        DO II = 1, IIMX
            CELLs(II)%flowVelocity(:) = 0.0
            num_node = size(CELLs(II)%nodeID)
            do n = 1, num_node
                ID = CELLs(II)%nodeID(n)
                CELLs(II)%flowVelocity(:) = CELLs(II)%flowVelocity(:) + pointVector(:,ID)
            end do
            CELLs(II)%flowVelocity(:) = CELLs(II)%flowVelocity(:) / num_node
        END DO
  
        ! do II = 1, size(celldata, dim=2)   !点データをセルデータに変換
              
        !     IF (celltype(II)==0) THEN
  
        !         celldata(:,II) = 0.25d0*(pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) &
        !             + pointdata(:, cell2node(3,II)) + pointdata(:, cell2node(4,II)))
  
        !     ELSE IF (celltype(II)==1) THEN
  
        !         celldata(:,II) = (pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) + pointdata(:, cell2node(3,II)) &
        !             + pointdata(:, cell2node(4,II)) + pointdata(:, cell2node(5,II)) + pointdata(:, cell2node(6,II))) / 6.0d0
  
        !     ELSE IF (celltype(II)==2) THEN
  
        !         celldata(:,II) = 0.20d0*(pointdata(:, cell2node(1,II)) + pointdata(:, cell2node(2,II)) &
        !             + pointdata(:, cell2node(3,II)) + pointdata(:, cell2node(4,II)) + pointdata(:, cell2node(5,II)))
  
        !     END IF
  
        ! end do
  
    end subroutine point2cellVelocity

    integer function get_mesh_info(name)
        character(*), intent(in) :: name
        
        select case(name)
            case('node')
                get_mesh_info = size(NODEs)

            case('cell')
                get_mesh_info = size(CELLs)

            case('tetra')
                get_mesh_info = count(CELLs(:)%typeName == 'tetra')

            case('prism')
                get_mesh_info = count(CELLs(:)%typeName == 'prism')

            case('pyramid')
                get_mesh_info = count(CELLs(:)%typeName == 'pyrmd')

            case default
                get_mesh_info = 0

        end select

    end function get_mesh_info
      
    subroutine deallocation_unstructuredGRID
        use vtkMesh_operator_m

        deallocate(CELLs)
        deallocate(NODEs)
        deallocate(BoundFACEs)
        if(FILE_TYPE=='VTK') call deallocation_VTK
    end subroutine deallocation_unstructuredGRID
    
end module unstructuredGrid_mod