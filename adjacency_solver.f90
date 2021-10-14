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

    ! INTEGER, allocatable :: NFN(:,:),NFC(:,:),NCF(:), NFNSUM(:)
    ! INTEGER, allocatable :: ICF(:,:)

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

    subroutine FACESET
        INTEGER II,JJJ, j, n, n_type, JJJMX
        integer, parameter :: num_halfFace_perCELL(3) = [4,5,5]
        integer, parameter :: IDtrans(4,5,3) = reshape([ &
                                1,2,3,0, 2,3,4,0, 3,4,1,0, 4,1,2,0, 0,0,0,0,&
                                1,2,3,0, 4,5,6,0, 1,2,4,5, 2,3,5,6, 3,1,6,4,&
                                5,1,2,0, 5,2,3,0, 5,3,4,0, 5,4,1,0, 1,2,3,4 ], shape(IDtrans))
                                                        

        num_halfFace = num_tetras*4 + num_prisms*5 + num_pyramids*5
        ! if(count(CELL_TYPE == 2) == 0)then
        !     JJTOTAL = 4*num_cells    !JJTOTAL:想定最大面数
        ! else
        !     JJTOTAL = 5*num_cells
        ! end if
    
        allocate(halfFACEs(num_halfFace))
        ! allocate(NFN(4,JJTOTAL), source = 0)
        ! allocate(NFC(2,JJTOTAL), source = -1)
        ! allocate(ICF(2,num_cells+1), source = 0)
        ! allocate(NFNSUM(JJTOTAL), source = 0)
            
            !     FACESET
        JJJ = 0   !単純面JJJ：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
        
        DO II = 1, num_cells
            select case(CELLs(II)%typeName)
                case('tetra')
                    n_type = 1
                case('prism')
                    n_type = 2
                case('pyrmd')
                    n_type = 3
                case default
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

        print*, 'ID_sum =', maxval(halfFACEs(:)%ID_sum)
            
    end subroutine  faceset
            !**************************************************************************************
        !**************************************************************************************
    subroutine facecheck
        INTEGER AA,BB,flag, DN,LL,LLX,FDN, JJJ,JJJ2,JJJa,JJJa2, numnode
        integer,allocatable :: JFS(:),JFL(:)!, sameface(:)
        !=======================================================================
        print*,'START-FACE CHECK!'
            !同一面の探索、真の面数算出
    
        FDN = num_halfFace/10000 + 1   !分割数目安（単純面数が多いほど分割数も多くなる）
    
        allocate(JFS(num_halfFace))
        allocate(JFL(2*FDN), source=0)
          
        DN = 3*num_nodes/FDN + 1   !分割幅（予想される最大節点番号和を分割数だけ分割）
        if(DN <= 0) then
            print*, 'ERROR_DN', DN
        end if
          
        JJJa = 1
        LL = 0
        do while(JJJa <= num_halfFace)  !面をグループに分けるループ（目的：JJJMXの2重ループを避けること）
            LL = LL +1
            JFL(LL) = JJJa
            do JJJ = 1, num_halfFace
                if((halfFACEs(JJJ)%ID_sum >= (LL-1)*DN).and.(halfFACEs(JJJ)%ID_sum < LL*DN)) then
                JFS(JJJa) = JJJ
                JJJa = JJJa + 1
                end if
            end do
        !print*,'LL,JFL=', LL, JFL(LL)
        end do
        LLX = LL
        print*, 'LLX=', LLX
        JFL(LLX+1) = JJJa
    
        if((JJJa-1) /= num_halfFace) then
            print*, 'JJJa_ERROR:', JJJa-1, num_halfFace
            stop
        end if
          
        ! allocate(sameface(JJJMX), source=0)
        num_BoundFaces = 0
        !$omp parallel do private(JJJ, JJJ2, flag, AA, BB, numnode) reduction(+:num_BoundFaces)
          
        LLloop : do LL = 1, LLX
        print*, 'CHECK_LL:', LL, '/', LLX
            face1 : do JJJa = JFL(LL), JFL(LL+1) -1
                JJJ = JFS(JJJa)
                if(halfFACEs(JJJ)%ownerID(2) >= 0) cycle face1 !面共有セル探索が済んでいる場合スキップ
                if(halfFACEs(JJJ)%nodeID(4) == 0) then
                    numnode = 3   !三角形面
                else
                    numnode = 4   !四角形面
                end if
        
                face2 :do JJJa2 = JJJa +1, JFL(LL+1) -1
                    JJJ2 = JFS(JJJa2)
                    if(halfFACEs(JJJ2)%ownerID(2) >= 0) cycle face2 !面共有セル探索が済んでいる場合スキップ
                    if((halfFACEs(JJJ2)%nodeID(4) > 0).and.(numnode==3)) cycle face2 !三角形面を注目中に四角形面が現れればスキップ
            
                    flag = 0
                    do AA = 1, numnode
                        do BB = 1, numnode
                            if(halfFACEs(JJJ)%nodeID(AA) == halfFACEs(JJJ2)%nodeID(BB)) flag = flag + 1  !点が一致すればカウント
                        end do
                    end do
            
                    if(flag < numnode) cycle face2 !節点数と同じ回数一致しなければスキップ（必要十分条件）
            
                    !ここまでくれば同一の2面発見
            
                    halfFACEs(JJJ)%ownerID(2) = halfFACEs(JJJ2)%ownerID(1)
                    halfFACEs(JJJ2)%ownerID(2) = halfFACEs(JJJ)%ownerID(1)
                    ! NFC(2,JJJ) = NFC(1,JJJ2)
                    ! NFC(2,JJJ2) = NFC(1,JJJ)
                    ! sameface(JJJ) = JJJ2
                    ! sameface(JJJ2) = JJJ 
            
                    cycle face1 !共有セルが見つかったので次の面へ
        
                end do face2
        
                halfFACEs(JJJ)%ownerID(2) = 0 !共有セルが見つからなかった（境界面）
                num_BoundFaces = num_BoundFaces + 1  !境界面カウント
        
            end do face1
    
        end do LLloop
          
        !$omp end parallel do
        
        ! allocate(NCF(JJJMX), source = 0)
        
        ! JJ = 0
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
        print*,'BC FACES',num_BoundFaces  
        print*,'END-FACE CHECK!'

    END subroutine facecheck

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

        call set_GRID_INFO  !要素数等の取得

        if(num_cells <= 0) then
            print*, 'ERROR_num_cells', num_cells
            stop
        end if

        call faceset    !面情報のセッティング
        call facecheck  !同一面のチェック

        call find_nodeID_onBoundFace   !境界面情報

        call solve_adjacency   !セル隣接情報

        call deallocation_adjacent  !配列解放

    end subroutine solve_adjacentInformation

END MODULE adjacent_information