MODULE adjacent_information
    use unstructuredGrid_mod
    IMPLICIT NONE
    integer num_nodes, num_cells, num_tetras, num_prisms, num_pyramids
    INTEGER JJMX, JJJMX, JJTOTAL
    integer num_squares, num_BoundFaces
    INTEGER, allocatable :: NFN(:,:),NFC(:,:),NCF(:), NFNSUM(:)
    INTEGER, allocatable :: ICF(:,:)

    contains

    subroutine set_GRID_INFO

        num_nodes = get_mesh_info('node')
        num_cells = get_mesh_info('cell')
        num_tetras = get_mesh_info('tetra')
        num_prisms = get_mesh_info('prism')
        num_pyramids = get_mesh_info('pyramid')

    end subroutine set_GRID_INFO

    subroutine FACESET
        INTEGER II,JJJ, face, node
        integer, parameter :: numface(3) = [4,5,5]
        integer, parameter :: ICNtrans(4,5,3) = reshape([ &
                                1,2,3,0, 2,3,4,0, 3,4,1,0, 4,1,2,0, 0,0,0,0,&
                                1,2,3,0, 4,5,6,0, 1,2,4,5, 2,3,5,6, 3,1,6,4,&
                                5,1,2,0, 5,2,3,0, 5,3,4,0, 5,4,1,0, 1,2,3,4 ], shape(ICNtrans))
                                                        
        num_squares = 0

        if(count(CELL_TYPE == 2) == 0)then
            JJTOTAL = 4*num_cells    !JJTOTAL:想定最大面数
        else
            JJTOTAL = 5*num_cells
        end if
    
        allocate(NFN(4,JJTOTAL), source = 0)
        allocate(NFC(2,JJTOTAL), source = -1)
        allocate(ICF(2,num_cells+1), source = 0)
        allocate(NFNSUM(JJTOTAL), source = 0)
            
            !     FACESET
        JJJ = 0   !単純面JJJ：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
        
        DO II = 1, num_cells
            
            do face = 1, numface(CELL_TYPE(II)+1)
                JJJ = JJJ + 1
                do node = 1, 3
                    NFN(node,JJJ) = ICN(ICNtrans(node, face, CELL_TYPE(II)+1), II)
                end do
                if(ICNtrans(4, face, CELL_TYPE(II)+1) > 0) then
                    NFN(4,JJJ) = ICN(ICNtrans(4, face, CELL_TYPE(II)+1), II) !四角形面なら4点目を代入
                    num_squares = num_squares + 1   !四角形面カウント
                end if
                NFC(1,JJJ) = II
                NFNSUM(JJJ) = NFN(1,JJJ) + NFN(2,JJJ) + NFN(3,JJJ) + NFN(4,JJJ)
                if (face == 1) then
                    ICF(1,II) = JJJ
                else if (face == numface(CELL_TYPE(II)+1)) then
                    ICF(2,II) = JJJ
                end if
            end do
            
        END DO    
              
        JJJMX = JJJ
        
        if (JJJMX == (num_tetras*4 + num_prisms*5 + num_pyramids*5)) then
            print*, 'JJJMX / JJTOTAL =', JJJMX, '/', JJTOTAL
        else
            print*, 'JJJMX_ERROR:', JJJMX, num_tetras, num_prisms, num_pyramids
            stop
        end if
    
        if(num_squares == (num_prisms*3 + num_pyramids)) then
            print*, 'num_squares=', num_squares
        else
            print*, 'SQUARE_ERROR:', num_squares, num_prisms, num_pyramids
            stop
        end if

        print*, 'NFNSUM=', maxval(NFNSUM)
            
    end subroutine  faceset
            !**************************************************************************************
        !**************************************************************************************
    subroutine facecheck
        INTEGER JJ,AA,BB,flag, DN,LL,LLX,FDN, JJJ,JJJ2,JJJa,JJJa2, numnode
        integer,allocatable :: JFS(:),JFL(:), sameface(:)
        !=======================================================================
        print*,'START-FACE CHECK!'
            !同一面の探索、真の面数算出
    
        FDN = JJJMX/10000 + 1   !分割数目安（単純面数が多いほど分割数も多くなる）
    
        allocate(JFS(JJJMX))
        allocate(JFL(2*FDN), source=0)
          
        DN = 3*num_nodes/FDN + 1   !分割幅（予想される最大節点番号和を分割数だけ分割）
        if(DN <= 0) then
            print*, 'ERROR_DN', DN
        end if
          
        JJJa = 1
        LL = 0
        do while(JJJa <= JJJMX)  !面をグループに分けるループ（目的：JJJMXの2重ループを避けること）
            LL = LL +1
            JFL(LL) = JJJa
            do JJJ = 1, JJJMX
                if((NFNSUM(JJJ)>=(LL-1)*DN).and.(NFNSUM(JJJ)<LL*DN)) then
                JFS(JJJa) = JJJ
                JJJa = JJJa + 1
                end if
            end do
        !print*,'LL,JFL=', LL, JFL(LL)
        end do
        LLX = LL
        print*, 'LLX=', LLX
        JFL(LLX+1) = JJJa
    
        if((JJJa-1) /= JJJMX) then
            print*, 'JJJa_ERROR:', JJJa-1, JJJMX
            stop
        end if
          
          
        allocate(sameface(JJJMX), source=0)
        num_BoundFaces = 0
        !$omp parallel do private(JJJ, JJJ2, flag, AA, BB, numnode) reduction(+:num_BoundFaces)
          
        LLloop : do LL = 1, LLX
        print*, 'CHECK_LL:', LL, '/', LLX
            face1 : do JJJa = JFL(LL), JFL(LL+1) -1
                JJJ = JFS(JJJa)
                if(NFC(2,JJJ) >= 0) cycle face1 !面共有セル探索が済んでいる場合スキップ
                if(NFN(4,JJJ) == 0) then
                    numnode = 3   !三角形面
                else
                    numnode = 4   !四角形面
                end if
        
                face2 :do JJJa2 = JJJa +1, JFL(LL+1) -1
                    JJJ2 = JFS(JJJa2)
                    if(NFC(2,JJJ2) >= 0) cycle face2 !面共有セル探索が済んでいる場合スキップ
                    if((NFN(4,JJJ2) > 0).and.(numnode==3)) cycle face2 !三角形面を注目中に四角形面が現れればスキップ
            
                    flag = 0
                    do AA = 1, numnode
                        do BB = 1, numnode
                            if(NFN(AA,JJJ) == NFN(BB,JJJ2)) flag = flag + 1  !点が一致すればカウント
                        end do
                    end do
            
                    if(flag < numnode) cycle face2 !節点数と同じ回数一致しなければスキップ（必要十分条件）
            
                    !ここまでくれば同一の2面発見
            
                    NFC(2,JJJ) = NFC(1,JJJ2)
                    NFC(2,JJJ2) = NFC(1,JJJ)
                    sameface(JJJ) = JJJ2
                    sameface(JJJ2) = JJJ 
            
                    cycle face1 !共有セルが見つかったので次の面へ
        
                end do face2
        
                NFC(2,JJJ) = 0 !共有セルが見つからなかった（境界面）
                num_BoundFaces = num_BoundFaces + 1  !境界面カウント
        
            end do face1
    
        end do LLloop
          
        !$omp end parallel do
        
        allocate(NCF(JJJMX), source = 0)
        
        JJ = 0
        DO JJJ = 1, JJJMX
            if ((NFC(1,JJJ) > num_cells).or.(NFC(2,JJJ) > num_cells).or.(NFC(1,JJJ) <= -1).or.(NFC(2,JJJ) <= -1)) then  !エラー条件
                print*, 'NFC_ERROR:', JJJ, NFC(1,JJJ), NFC(2,JJJ), '/', num_cells
                stop
            
            else if ((NFC(2,JJJ) > NFC(1,JJJ)) .OR. (NFC(2,JJJ) == 0)) THEN !この条件により同一の2面のうちの一方が棄却される
                JJ = JJ + 1 !真の面数カウント
                NFN(1,JJ) = NFN(1,JJJ)  !常に JJ<=JJJ であり、棄却された面だけ前に詰めて代入
                NFN(2,JJ) = NFN(2,JJJ)
                NFN(3,JJ) = NFN(3,JJJ)
                NFN(4,JJ) = NFN(4,JJJ)
            
                NFC(1,JJ) = NFC(1,JJJ)
                NFC(2,JJ) = NFC(2,JJJ)
            
                NCF(JJJ) = JJ
                if(sameface(JJJ) > 0) NCF(sameface(JJJ)) = JJ   !同一面が存在すればその面に対しても処理
            
            end if
        END DO
        JJMX = JJ  !これが真の面数
          
        print*,'true_FACES',JJMX
        print*,'BC FACES',num_BoundFaces  
        print*,'END-FACE CHECK!'
          !
        
          
    END subroutine facecheck

    subroutine find_nodeID_onBoundFace
        INTEGER II,JJ,JB
        integer, allocatable :: nodeID_onBoundFace(:,:)
  
        allocate(NoB(num_cells), source=0)
        allocate(ICB(4,num_cells), source=0) !II consists Boundaries

        allocate(nodeID_onBoundFace(3, num_BoundFaces))
        JB = 0

        DO JJ = 1, JJMX
            IF (NFC(2,JJ) /= 0) cycle   !境界面以外はスルー
            JB = JB +1
            II = NFC(1,JJ) !JJが属する要素番号
            NoB(II) = NoB(II) +1 !IIが所有する境界面数カウント
            ICB(NoB(II),II) = JB !IIが所有する境界面番号
    
            nodeID_onBoundFace(1:3, JB) = NFN(1:3,JJ)

        END DO

        call set_nodeID_onBoundFace(nodeID_onBoundFace)
       
    END subroutine find_nodeID_onBoundFace
        
    subroutine solve_adjacency
        integer II, JJ, JJJ, numnext, cnt
        integer, parameter :: LF = 1
        integer, allocatable :: num_adjacent(:), adjacentCellID(:,:)
        
        allocate(num_adjacent(num_cells))

        if(num_prisms==0)then
            allocate(adjacentCellID(4, num_cells))
        else
            allocate(adjacentCellID(5, num_cells))
        end if
    
        do II = 1, num_cells
            if(CELL_TYPE(II) == 0) then !テトラ
                numnext = 4 !隣接セル数
            else if(CELL_TYPE(II) == 1) then
                numnext = 5
            else if(CELL_TYPE(II) == 2) then
                numnext = 5
            else
                print*, 'CELL_TYPE error:', CELL_TYPE(II)
                stop
            end if
            num_adjacent(II) = numnext
    
            cnt = 1
            do JJJ = ICF(1,II), ICF(1,II) + numnext -1
                JJ = NCF(JJJ)
                if((JJ <= 0).or.(JJ > JJMX)) then
                    print*, 'NCF_ERROR', JJ, JJJ, II
                    stop
                end if
        
                if (NFC(1,JJ) == II) then
                    adjacentCellID(cnt,II) = NFC(2,JJ)
                else
                    adjacentCellID(cnt,II) = NFC(1,JJ)
                end if
                cnt = cnt + 1
            end do
    
        end do

        call set_adjacency(num_adjacent, adjacentCellID)

    end subroutine solve_adjacency

    subroutine deallocation_adjacent

        print*, 'Deallocation AdjacentInfo'

        if(allocated(NFN)) deallocate(NFN)
        if(allocated(NFC)) deallocate(NFC)
        if(allocated(NCF)) deallocate(NCF)
        if(allocated(NFNSUM)) deallocate(NFNSUM)
        if(allocated(ICF)) deallocate(ICF)

    end subroutine deallocation_adjacent

    subroutine solve_adjacentInformation

        call set_GRID_INFO  !要素番号等の取得

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