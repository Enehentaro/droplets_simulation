MODULE adjacent_information
    use unstructuredGrid_mod
    IMPLICIT NONE
    integer num_nodes, num_cells, num_tetras, num_prisms, num_pyramids
    INTEGER num_halfFace
    integer num_BoundFaces

    type halfFace_t
        integer :: nodeID(4) = 0, pairID = 0, ID_sum = 0, ownerID(2) = -1
    end type halfFace_t

    type(halfFace_t), allocatable :: halfFACEs(:)

    contains

    subroutine set_GRID_INFO

        num_nodes = get_mesh_info('node')
        num_cells = get_mesh_info('cell')
        num_tetras = get_mesh_info('tetra')
        num_prisms = get_mesh_info('prism')
        num_pyramids = get_mesh_info('pyramid')

        print*, 'NODE :', num_nodes
        print*, 'CELL :', num_cells
        print*, 'tet,pri,pyr :', num_tetras, num_prisms, num_pyramids

    end subroutine set_GRID_INFO

    subroutine set_halfFACEs
        INTEGER II,JJJ, j, n, n_type, JJJMX
        integer, parameter :: num_halfFace_perCELL(3) = [4,5,5]
        integer, parameter :: IDtrans(4,5,3) = reshape([ &
                                1,2,3,0, 2,3,4,0, 3,4,1,0, 4,1,2,0, 0,0,0,0,&
                                1,2,3,0, 4,5,6,0, 1,2,4,5, 2,3,5,6, 3,1,6,4,&
                                5,1,2,0, 5,2,3,0, 5,3,4,0, 5,4,1,0, 1,2,3,4 ], shape(IDtrans))
                                                        

        num_halfFace = num_tetras*4 + num_prisms*5 + num_pyramids*5   !半面数：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
    
        allocate(halfFACEs(num_halfFace))
            
        JJJ = 0
        
        DO II = 1, num_cells
            select case(CELLs(II)%typeName)
                case('tetra')
                    n_type = 1
                case('prism')
                    n_type = 2
                case('pyrmd')
                    n_type = 3
                case default
                    n_type = -1
            end select

            do j = 1, num_halfFace_perCELL(n_type)
                JJJ = JJJ + 1
                do n = 1, 3
                    ! print*, JJJ, n, j, n_type, II
                    halfFACEs(JJJ)%nodeID(n) = CELLs(II)%nodeID(IDtrans(n, j, n_type))
                    ! NFN(node,JJJ) = ICN(ICNtrans(node, face, CELL_TYPE(II)+1), II)
                end do
                if(IDtrans(4, j, n_type) > 0) then
                    halfFACEs(JJJ)%nodeID(4) = CELLs(II)%nodeID(IDtrans(4, j, n_type))
                    ! NFN(4,JJJ) = ICN(ICNtrans(4, face, CELL_TYPE(II)+1), II) !四角形面なら4点目を代入
                end if
                halfFACEs(JJJ)%ownerID(1) = II
                halfFACEs(JJJ)%ID_sum = sum(halfFACEs(JJJ)%nodeID(:))
                ! NFNSUM(JJJ) = NFN(1,JJJ) + NFN(2,JJJ) + NFN(3,JJJ) + NFN(4,JJJ)
            end do
            
        END DO    
              
        JJJMX = JJJ
        
        if (JJJMX == (num_tetras*4 + num_prisms*5 + num_pyramids*5)) then
            print*, 'JJJMX / num_halfFace =', JJJMX, '/', num_halfFace
        else
            print*, 'JJJMX_ERROR:', JJJMX, num_tetras, num_prisms, num_pyramids
            stop
        end if
            
    end subroutine  set_halfFACEs
    
    subroutine check_FACEs
        use terminalControler_m
        INTEGER match, width, numNode, faceID, groupID,num_group, maxID_sum
        integer checkCounter, faceCounter, i,j, faceID1,faceID2, num_face, k,l
        integer, allocatable :: faceID_array(:)
        type faceGroup_t
            integer, allocatable :: faceID(:)
        end type faceGroup_t
        type(faceGroup_t), allocatable :: faceGroup(:)
        ! real time1, time2

        ! call cpu_time(time1)

        print*,'START-FACE CHECK!' !同一面の探索
            
        num_group = num_halfFace/5000 + 1   !面グループ数（1グループ数に約5000枚面が入るようにする）（この値は経験則）
    
        allocate(faceGroup(num_group))
        allocate(faceID_array(num_halfFace), source=0)

        maxID_sum = maxval(halfFACEs(:)%ID_sum)
        width = maxID_sum/num_group + 1   !1グループの幅（最大節点番号和を面グループ数で割る）
          
        checkCounter = 0
        call set_formatTC('("DIVIDE halfFace [ #group : ",i6," / ",i6," ]")')
        do groupID = 1, num_group  !面をグループに分ける
            call print_sameLine([groupID, num_group])
            faceID_array(:) = 0
            faceCounter = 1
            do faceID = 1, num_halfFace
                if((halfFACEs(faceID)%ID_sum > (groupID-1)*width).and.(halfFACEs(faceID)%ID_sum <= groupID*width)) then
                    faceID_array(faceCounter) = faceID
                    faceCounter = faceCounter + 1
                    checkCounter = checkCounter + 1
                end if
            end do
            faceGroup(groupID)%faceID = faceID_array(1 : faceCounter-1)
        end do

        if(checkCounter /= num_halfFace) then
            print*, 'faceCounter_ERROR:', checkCounter, num_halfFace
            stop
        end if
          
        num_BoundFaces = 0
        call set_formatTC('("CHECK halfFace [ #group : ",i6," / ",i6," ]")')
        !$omp parallel do private(faceID1,faceID2, match, k,l, num_face,numNode) reduction(+:num_BoundFaces)
        do groupID = 1, num_group
            call print_sameLine([groupID, num_group])
            num_face = size(faceGroup(groupID)%faceID)

            face1 : do i = 1, num_face
                faceID1 = faceGroup(groupID)%faceID(i)
                if(halfFACEs(faceID1)%ownerID(2) >= 0) cycle face1 !面共有セル探索が済んでいる場合スキップ
                if(halfFACEs(faceID1)%nodeID(4) <= 0) then
                    numNode = 3   !三角形面
                else
                    numNode = 4   !四角形面
                end if
        
                face2 :do j = i+1, num_face
                    faceID2 = faceGroup(groupID)%faceID(j)
                    if(halfFACEs(faceID2)%ownerID(2) >= 0) cycle face2 !面共有セル探索が済んでいる場合スキップ
                    if((numNode==3).and.(halfFACEs(faceID2)%nodeID(4) > 0)) cycle face2 !三角形面を注目中に四角形面が現れればスキップ
            
                    match = 0
                    do k = 1, numNode
                        do l = 1, numNode
                            if(halfFACEs(faceID1)%nodeID(k) == halfFACEs(faceID2)%nodeID(l)) match = match + 1  !点IDが一致すればカウント
                        end do
                    end do
            
                    if(match < numNode) cycle face2 !節点数と同じ回数一致しなければスキップ（必要十分条件）
            
                    !ここまでくれば同一の2面発見
            
                    halfFACEs(faceID1)%ownerID(2) = halfFACEs(faceID2)%ownerID(1)
                    halfFACEs(faceID2)%ownerID(2) = halfFACEs(faceID1)%ownerID(1)
            
                    cycle face1 !共有セルが見つかったので次の面へ
        
                end do face2
        
                halfFACEs(faceID1)%ownerID(2) = 0 !共有セルが見つからなかった:境界面
                num_BoundFaces = num_BoundFaces + 1  !境界面カウント
        
            end do face1
    
        end do
        !$omp end parallel do
        
        ! allocate(NCF(JJJMX), source = 0)
        
        ! JJ = 0    !真の面数算出
        ! DO JJJ = 1, JJJMX
        !     if ((NFC(1,JJJ) > num_cells).or.(NFC(2,JJJ) > num_cells).or.(NFC(1,JJJ) <= -1).or.(NFC(2,JJJ) <= -1)) then  !エラー条件
        !         print*, 'NFC_ERROR:', JJJ, NFC(1,JJJ), NFC(2,JJJ), '/', num_cells
        !         stop
            
        !     else if ((NFC(2,JJJ) > NFC(1,JJJ)) .OR. (NFC(2,JJJ) == 0)) THEN !この条件により同一の2面のうちの一方が棄却される
        !         JJ = JJ + 1 !真の面数カウント
        !         NFN(1,JJ) = NFN(1,JJJ)  !常に JJ<=JJJ であり、棄却された面だけ前に詰めて代入
        !         NFN(2,JJ) = NFN(2,JJJ)
        !         NFN(3,JJ) = NFN(3,JJJ)
        !         NFN(4,JJ) = NFN(4,JJJ)
            
        !         NFC(1,JJ) = NFC(1,JJJ)
        !         NFC(2,JJ) = NFC(2,JJJ)
            
        !         NCF(JJJ) = JJ
        !         if(sameface(JJJ) > 0) NCF(sameface(JJJ)) = JJ   !同一面が存在すればその面に対しても処理
            
        !     end if
        ! END DO
        ! JJMX = JJ  !これが真の面数
          
        ! print*,'true_FACES',JJMX
        print*,'# Boundary Face =', num_BoundFaces
        print*,'END-FACE CHECK!'

        ! call cpu_time(time2)
        ! print*, time2 - time1
        ! stop

    END subroutine check_FACEs

    subroutine find_nodeID_onBoundFace
        INTEGER II,JJJ,JB
        integer, allocatable :: NoB(:), ICB(:,:)
  
        allocate(NoB(num_cells), source=0)
        allocate(ICB(5, num_cells), source=0) !II consists Boundaries

        allocate(BoundFACEs(num_BoundFaces))
        JB = 0

        DO JJJ = 1, num_halfFace
            IF (halfFACEs(JJJ)%ownerID(2) /= 0) cycle   !境界面以外はスルー
            JB = JB +1
            II = halfFACEs(JJJ)%ownerID(1) !JJが属する要素番号
            NoB(II) = NoB(II) +1 !IIが所有する境界面数カウント
            ICB(NoB(II),II) = JB !IIが所有する境界面番号
    
            BoundFACEs(JB)%nodeID(1:3) = halfFACEs(JJJ)%nodeID(1:3)

        END DO

        do II = 1, num_cells
            CELLs(II)%boundFaceID = ICB(1:NoB(II), II)
        end do

        ! call set_nodeID_onBoundFace(nodeID_onBoundFace)
       
    END subroutine find_nodeID_onBoundFace
        
    subroutine solve_adjacency
        integer II, JJJ!, JJ, numnext, cnt
        integer, parameter :: LF = 1
        integer, allocatable :: num_adjacent(:), adjacentCellID(:,:)
        
        allocate(num_adjacent(num_cells), source=0)

        if(num_prisms==0)then
            allocate(adjacentCellID(4, num_cells))
        else
            allocate(adjacentCellID(5, num_cells))
        end if

        do JJJ = 1, num_halfFace
            if (halfFACEs(JJJ)%ownerID(2) <= 0) cycle !境界面はスルー

            II = halfFACEs(JJJ)%ownerID(1)
            num_adjacent(II) = num_adjacent(II) + 1
            adjacentCellID(num_adjacent(II), II) = halfFACEs(JJJ)%ownerID(2)

        end do

        do II = 1, num_cells
            CELLs(II)%adjacentCellID = adjacentCellID(1:num_adjacent(II), II)
        end do

    
        ! do II = 1, num_cells
        !     if(CELL_TYPE(II) == 0) then !テトラ
        !         numnext = 4 !隣接セル数
        !     else if(CELL_TYPE(II) == 1) then
        !         numnext = 5
        !     else if(CELL_TYPE(II) == 2) then
        !         numnext = 5
        !     else
        !         print*, 'CELL_TYPE error:', CELL_TYPE(II)
        !         stop
        !     end if
        !     num_adjacent(II) = numnext
    
        !     cnt = 1
        !     do JJJ = ICF(1,II), ICF(1,II) + numnext -1
        !         JJ = NCF(JJJ)
        !         if((JJ <= 0).or.(JJ > JJMX)) then
        !             print*, 'NCF_ERROR', JJ, JJJ, II
        !             stop
        !         end if
        
        !         if (NFC(1,JJ) == II) then
        !             adjacentCellID(cnt,II) = NFC(2,JJ)
        !         else
        !             adjacentCellID(cnt,II) = NFC(1,JJ)
        !         end if
        !         cnt = cnt + 1
        !     end do
    
        ! end do

        ! call set_adjacency(num_adjacent, adjacentCellID)

    end subroutine solve_adjacency

    subroutine deallocation_adjacent

        print*, 'Deallocation AdjacentInfo'

        deallocate(halfFACEs)

    end subroutine deallocation_adjacent

    subroutine solve_adjacentInformation

        print*, '~ SolvingAdjacency is Required. ~'

        call set_GRID_INFO  !要素数等の取得

        if(num_cells <= 0) then
            print*, 'ERROR_num_cells', num_cells
            stop
        end if

        call set_halfFACEs    !面情報のセッティング
        call check_FACEs  !同一面のチェック

        call find_nodeID_onBoundFace   !境界面情報

        call solve_adjacency   !セル隣接情報

        call deallocation_adjacent  !配列解放

    end subroutine solve_adjacentInformation

END MODULE adjacent_information