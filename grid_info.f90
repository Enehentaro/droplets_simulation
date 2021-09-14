!==============================================================================================================
MODULE grid_information
    IMPLICIT NONE
    INTEGER JJMX,KKMX,JJJMX,LBX, JBMX
    INTEGER IITPRMX,JJTPRMX
    integer :: IIMX,IITETMX,IIPRSMX,IIPYRMX, INP=0, numfiles=1, SQUARES=0, digits=0
    INTEGER JJTOTAL
    INTEGER, allocatable :: NFN(:,:),NFC(:,:),NCF(:), NFNSUM(:)
    INTEGER, allocatable :: ICF(:,:),JUDTP(:)
    INTEGER, allocatable :: ICN(:,:)
    INTEGER, allocatable :: cell_to_cell(:,:), NoB(:), nearmax(:), ICB(:,:)
    DOUBLE PRECISION, allocatable :: R(:,:), UVW(:,:)

    contains

!**************************************************************************************
    subroutine readfile(FNAME)
        INTEGER II,II2,KK, IIMX2, n_unit, L
        CHARACTER(*), intent(in) :: FNAME
        CHARACTER*15 AAA
        INTEGER num
        !=======================================================================
        
        !-------データ読込み---------------------------------------------------------------------------
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,*)AAA,KKMX
            print*, 'KKMX=', KKMX
            allocate(R(3,KKMX))
            DO KK = 1, KKMX
                read(n_unit,*)R(1,KK),R(2,KK),R(3,KK)
            END DO
            read(n_unit,'()')
            read(n_unit,*)AAA,IIMX2
            print*, 'IIMX2=', IIMX2
        
            allocate(ICN(6,IIMX2), source = -1)
            allocate(JUDTP(IIMX2), source = 0)
        
            II = 0
    
            DO II2 = 1, IIMX2
                read(n_unit, fmt='(I12)', advance='no')num !ここのフォーマットは、読み込むファイルに合わせて変更の必要があるかも
                
                if((num>=4).and.(num<=6)) then
                    II = II + 1
                    read(n_unit,*) (ICN(L,II), L=1,num)

                    select case(num)
                        case(4)
                            JUDTP(II) = 0
                            IITETMX=IITETMX+1   
                        case(6)
                            JUDTP(II) = 1
                            IIPRSMX=IIPRSMX+1
                        case(5)
                            JUDTP(II) = 2
                            IIPYRMX=IIPYRMX+1
                    end select
        
                else
                    read(n_unit,'()')  !面要素はスルー
        
                end if
    
            END DO
        
        close(n_unit)
        
        ICN(:,:) = ICN(:,:) + 1 !VTKを読んだ後は節点ID+1

        if(II == (IITETMX+IIPRSMX+IIPYRMX)) then
            IIMX = II
        else
            print*, 'IIMX_ERROR:', IIMX,IITETMX+IIPRSMX+IIPYRMX
            stop
        end if
        print*, 'KKMX=',KKMX
        print*, 'IIMX=',IIMX
        print*, 'II consists of',IITETMX,IIPRSMX,IIPYRMX
        
        IITPRMX=IITETMX+IIPRSMX
             
    end subroutine readfile
        !**************************************************************************************
        
        !**************************************************************************************
    subroutine readINP(FNAME)
        CHARACTER(*), intent(in) :: FNAME
        INTEGER II,II2,KK,AAmax, IIMX2, AA, n_unit
        character(6) cellshape
        double precision, allocatable :: UVWK(:,:)
            !=======================================================================
        IITETMX = 0
        IIPRSMX = 0
        IIPYRMX = 0
        print*, 'READINP:', FNAME
        
        open(newunit=n_unit,FILE=FNAME,STATUS='OLD')
            read(n_unit,*)KKMX,IIMX2
            print*, 'KKMX,IIMX2=',KKMX,IIMX2
        
            if(.not.allocated(R)) allocate(R(3,KKMX))
            if(.not.allocated(ICN)) allocate(ICN(6,IIMX2))
            if(.not.allocated(JUDTP)) allocate(JUDTP(IIMX2))
        
            allocate(UVWK(KKMX,3))
               
            DO KK = 1, KKMX
                read(n_unit,*)AA,R(1,KK),R(2,KK),R(3,KK)
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
                    JUDTP(II) = 0
                    IITETMX = IITETMX +1
                    read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II)    
                    if (ICN(1,II)==0.or.ICN(4,II)==0) print*, 'ICN_WARNING_tet:', ICN(:,II)
                    
                    ELSE IF(cellshape=='prism') THEN
                    JUDTP(II) = 1
                    IIPRSMX = IIPRSMX +1
                    read(n_unit,*)ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II),ICN(5,II),ICN(6,II)
                    if (ICN(1,II)==0.or.ICN(6,II)==0) print*, 'ICN_WARNING_prism:', ICN(:,II)
    
                    ELSE IF(cellshape=='pyr') THEN
                    JUDTP(II) = 2
                    IIPYRMX = IIPYRMX +1
                    read(n_unit,*)ICN(5,II),ICN(1,II),ICN(2,II),ICN(3,II),ICN(4,II) !INPは最初が山頂点であり、VTKでは最後が山頂点のため、読み込む順がこうなる。
                    if (ICN(1,II)==0.or.ICN(5,II)==0) print*, 'ICN_WARNING_pyr:', ICN(:,II)
    
                    end if
    
                ELSE
                    read(n_unit,'()')  !テトラでもプリズムでもピラミッドでもないならスルー
        
                ENDIF
        
            END DO
        
            if (II==(IITETMX + IIPRSMX + IIPYRMX)) then
            IIMX = II
            print*, 'IIMX=', IIMX
            print*, 'Tetra,Prism,Pyramid=', IITETMX, IIPRSMX, IIPYRMX
            else
            print*, 'IIMX_ERROR', II, (IITETMX + IIPRSMX + IIPYRMX)
            stop
            end if
            
            if(.not.allocated(UVW)) allocate(UVW(IIMX,3))
        
            read(n_unit,*)AAmax
            read(n_unit,'()')
            DO II = 1,AAmax
            read(n_unit,'()')
            END DO
            DO KK = 1, KKMX
            read(n_unit,*)AA,UVWK(KK,1),UVWK(KK,2),UVWK(KK,3)
            END DO
            
        close(n_unit)
        
        
        do II = 1, IIMX   !点データをセルデータに変換
        
            IF (JUDTP(II)==0) THEN

                UVW(II,1) = 0.25d0*(UVWK(ICN(1,II),1)+UVWK(ICN(2,II),1)+UVWK(ICN(3,II),1)+UVWK(ICN(4,II),1))
                UVW(II,2) = 0.25d0*(UVWK(ICN(1,II),2)+UVWK(ICN(2,II),2)+UVWK(ICN(3,II),2)+UVWK(ICN(4,II),2))
                UVW(II,3) = 0.25d0*(UVWK(ICN(1,II),3)+UVWK(ICN(2,II),3)+UVWK(ICN(3,II),3)+UVWK(ICN(4,II),3))
            ELSE IF (JUDTP(II)==1) THEN

                UVW(II,1) = (UVWK(ICN(1,II),1)+UVWK(ICN(2,II),1)+UVWK(ICN(3,II),1)+UVWK(ICN(4,II),1)&
                            +UVWK(ICN(5,II),1)+UVWK(ICN(6,II),1))/6.0d0
                UVW(II,2) = (UVWK(ICN(1,II),2)+UVWK(ICN(2,II),2)+UVWK(ICN(3,II),2)+UVWK(ICN(4,II),2)&
                            +UVWK(ICN(5,II),2)+UVWK(ICN(6,II),2))/6.0d0
                UVW(II,3) = (UVWK(ICN(1,II),3)+UVWK(ICN(2,II),3)+UVWK(ICN(3,II),3)+UVWK(ICN(4,II),3)&
                            +UVWK(ICN(5,II),3)+UVWK(ICN(6,II),3))/6.0d0
            ELSE IF (JUDTP(II)==2) THEN

                UVW(II,1) = 0.20d0*(UVWK(ICN(1,II),1)+UVWK(ICN(2,II),1)+UVWK(ICN(3,II),1)+UVWK(ICN(4,II),1)+UVWK(ICN(5,II),1))
                UVW(II,2) = 0.20d0*(UVWK(ICN(1,II),2)+UVWK(ICN(2,II),2)+UVWK(ICN(3,II),2)+UVWK(ICN(4,II),2)+UVWK(ICN(5,II),2))
                UVW(II,3) = 0.20d0*(UVWK(ICN(1,II),3)+UVWK(ICN(2,II),3)+UVWK(ICN(3,II),3)+UVWK(ICN(4,II),3)+UVWK(ICN(5,II),3))
            END IF
        
            !print*, 'UVW=', II, UVW(II,1),UVW(II,2),UVW(II,3)
        
        end do
        
        !deallocate(UVWK)
        
    end subroutine readINP
        !**************************************************************************************
        !**************************************************************************************
    subroutine FACESET
        INTEGER II,JJJ, face, node
        integer, parameter :: numface(3) = [4,5,5]
        integer, parameter :: ICNtrans(3,5,4) = reshape([1,2,3,0, 2,3,4,0, 3,4,1,0, 4,1,2,0, 0,0,0,0,&
                                                        1,2,3,0, 4,5,6,0, 1,2,4,5, 2,3,5,6, 3,1,6,4,&
                                                        5,1,2,0, 5,2,3,0, 5,3,4,0, 5,4,1,0, 1,2,3,4], [3,5,4], order=[3,2,1])
            
        print*, ICNtrans(1,1,:)
        print*, ICNtrans(2,1,:)
        print*, ICNtrans(3,1,:)
        
        if(IIPRSMX==0)then
            JJTOTAL = 4*IIMX    !JJTOTAL:想定最大面数
        else
            JJTOTAL = 5*IIMX
        end if
    
        allocate(NFN(4,JJTOTAL), source = 0)
        allocate(NFC(2,JJTOTAL), source = -1)
        allocate(ICF(2,IIMX+1), source = 0)
        allocate(NFNSUM(JJTOTAL), source = 0)
            
            !     FACESET
        JJJ = 0   !単純面JJJ：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
            
        DO II = 1, IIMX
            
            do face = 1, numface(JUDTP(II)+1)
                JJJ = JJJ + 1
                do node = 1, 3
                    NFN(node,JJJ) = ICN(ICNtrans(JUDTP(II)+1, face, node), II)
                end do
                if(ICNtrans(JUDTP(II)+1, face, 4) > 0) then
                    NFN(4,JJJ) = ICN(ICNtrans(JUDTP(II)+1, face, 4), II) !四角形面なら4点目を代入
                    SQUARES = SQUARES + 1   !四角形面カウント
                end if
                NFC(1,JJJ) = II
                NFNSUM(JJJ) = NFN(1,JJJ) + NFN(2,JJJ) + NFN(3,JJJ) + NFN(4,JJJ)
                if (face == 1) then
                    ICF(1,II) = JJJ
                else if (face == numface(JUDTP(II)+1)) then
                    ICF(2,II) = JJJ
                end if
            end do
            
        END DO    
              
        JJJMX = JJJ
        
        if (JJJMX == (IITETMX*4 + IIPRSMX*5 + IIPYRMX*5)) then
            print*, 'JJJMX / JJTOTAL =', JJJMX, '/', JJTOTAL
        else
            print*, 'JJJMX_ERROR:', JJJMX, IITETMX, IIPRSMX, IIPYRMX
            stop
        end if
    
        if(SQUARES == (IIPRSMX*3 + IIPYRMX)) then
            print*, 'SQUARES=', SQUARES
        else
            print*, 'SQUARE_ERROR:', SQUARES, IIPRSMX, IIPYRMX
            stop
        end if
            
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
          
        DN = 3*KKMX/FDN + 1   !分割幅（予想される最大節点番号和を分割数だけ分割）
          
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
        JBMX = 0
        !$omp parallel do private(JJJ, JJJ2, flag, AA, BB, numnode) reduction(+:JBMX)
          
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
                JBMX = JBMX + 1  !境界面カウント
        
            end do face1
    
        end do LLloop
          
        !$omp end parallel do
        
        allocate(NCF(JJJMX), source = 0)
        
        JJ = 0
        DO JJJ = 1, JJJMX
            if ((NFC(1,JJJ) > IIMX).or.(NFC(2,JJJ) > IIMX).or.(NFC(1,JJJ) <= -1).or.(NFC(2,JJJ) <= -1)) then  !エラー条件
                print*, 'NFC_ERROR:', JJJ, NFC(1,JJJ), NFC(2,JJJ), '/', IIMX
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
                JJJ2 = sameface(JJJ)    !JJJと同一の面番号
                NCF(JJJ2) = JJ
            
            end if
        END DO
        JJMX = JJ  !これが真の面数
          
        print*,'true_FACES',JJMX
        print*,'BC FACES',JBMX  
        print*,'END-FACE CHECK!'
          !
        
          
    END subroutine facecheck
          !**************************************************************************************
          !**************************************************************************************
           
        !**************************************************************************************
    subroutine boundaryset(dir)
        INTEGER II,JJ,JB, n_unit
        integer, parameter :: LF=1
        character(*), intent(in) :: dir
        character(99) FNAME
        !=======================================================================
        !-----BC SET-------------------------------------------------------
        allocate(NoB(IIMX), source=0)
        allocate(ICB(4,IIMX), source=0) !II consists Boundaries
        
        FNAME = trim(dir)//'/boundaries.txt'
        open(newunit=n_unit,FILE= FNAME, STATUS='REPLACE')
            write(n_unit,*) JBMX
    
            JB = 0
    
            DO JJ = 1, JJMX
                IF (NFC(2,JJ) /= 0) cycle   !境界面以外はスルー
                JB = JB +1
                II = NFC(1,JJ) !JJが属する要素番号
                NoB(II) = NoB(II) +1 !IIが所有する境界面数カウント
                ICB(NoB(II),II) = JB !IIが所有する境界面番号
        
                write(n_unit,*) NFN(1,JJ), NFN(2,JJ), NFN(3,JJ)
            END DO
    
        close(n_unit)
        
        if (JB /= JBMX) then
            print*, 'ERROR_JBMX', JB, JBMX
            stop
        end if
    
        print*, 'WRITEOUT:', FNAME
        
        !=======================================================================
    END subroutine boundaryset
        !**************************************************************************************
        
        !=======================================================================
    subroutine nextcell(dir)
        integer II, JJ, JJJ, numnext, JB, n_unit
        integer, parameter :: LF = 1
        character(*), intent(in) :: dir
        character(99) FNAME
    
        FNAME = trim(dir)//'/nextcell.txt'
          
        open(newunit=n_unit,FILE= FNAME, STATUS='REPLACE')
        
            write(n_unit,*) IIMX
        
            if(IIPRSMX==0)then
                write(n_unit,*) 4
            else
                write(n_unit,*) 5
            end if
        
            do II = 1, IIMX
                if(JUDTP(II) == 0) then !テトラ
                    numnext = 4 !隣接セル数
                else if(JUDTP(II) == 1) then
                    numnext = 5
                else if(JUDTP(II) == 2) then
                    numnext = 5
                end if
                write(n_unit, fmt='(I5)', advance='no') numnext
        
            do JJJ = ICF(1,II), ICF(1,II) + numnext -1
                JJ = NCF(JJJ)
                if((JJ <= 0).or.(JJ > JJMX)) then
                    print*, 'NCF_ERROR', JJ, JJJ, II
                    close(n_unit)
                    stop
                end if
        
                if (NFC(1,JJ) == II) then
                    write(n_unit, fmt='(I12)', advance='no') NFC(2,JJ)
                else
                    write(n_unit, fmt='(I12)', advance='no') NFC(1,JJ)
                end if
                
            end do
        
            write(n_unit,'()')  !改行
        
            end do
        
            DO II = 1, IIMX
                write(n_unit, fmt='(I4)', advance='no') NoB(II)  ! Number of Boundary
                do JB = 1, NoB(II)
                    write(n_unit, fmt='(I10)', advance='no') ICB(JB,II)
                end do
                write(n_unit,'()')  !改行
            END DO
        
        close(n_unit)
          
          
        print*, 'WRITEOUT:', FNAME
    end subroutine nextcell
          !=======================================================================
        
        !**************************************************************************************
    subroutine writeout
        INTEGER IITOTAL,LF, n_unit
        INTEGER II,KK
        character(20) FNAME
        !=======================================================================
        !-----WRITE OUT------
        IITOTAL=IITETMX*5+IIPRSMX*7+IIPYRMX*6
        
        FNAME = 'IGRID.vtk'
        
        LF = 1
        open(newunit=n_unit,FILE=FNAME,STATUS='REPLACE')
            write(n_unit,'(A)') '# vtk DataFile Version 2.0'
            write(n_unit,'(A)') 'FOR TEST'
            write(n_unit,'(A)') 'ASCII'
            write(n_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(A,I12,A)') 'POINTS ',KKMX,' float'
            DO KK = 1, KKMX
                write(n_unit,*)R(1,KK),R(2,KK),R(3,KK)
            END DO
            write(n_unit,'(A)')''
            write(n_unit,'(A,I12,2X,I12)') 'CELLS ',IIMX, IITOTAL
        
            DO II = 1, IIMX
                IF(JUDTP(II)==0)THEN
                    write(n_unit,'(5(I12,2X))')4,ICN(1,II)-1,ICN(2,II)-1,ICN(3,II)-1,ICN(4,II)-1
                ELSE IF(JUDTP(II)==1)THEN
                    write(n_unit,'(7(I12,2X))')6,ICN(1,II)-1,ICN(2,II)-1,ICN(3,II)-1,ICN(4,II)-1,ICN(5,II)-1,ICN(6,II)-1
                ELSE IF(JUDTP(II)==2)THEN
                    write(n_unit,'(6(I12,2X))')5,ICN(1,II)-1,ICN(2,II)-1,ICN(3,II)-1,ICN(4,II)-1,ICN(5,II)-1
                END IF
            END DO
        
            write(n_unit,'(A)')''
            write(n_unit,'(A,I12)') 'CELL_TYPES',IIMX
            DO II = 1, IIMX
                IF(JUDTP(II)==0)THEN
                    write(n_unit,'(I12)')10
                ELSE IF(JUDTP(II)==1)THEN
                    write(n_unit,'(I12)')13
                ELSE IF(JUDTP(II)==2)THEN
                    write(n_unit,'(I12)')14  
                END IF
            END DO
            write(n_unit,'(A)')''
        
            write(n_unit,'(A,I12)') 'CELL_DATA ',IIMX
            write(n_unit,'(A)') 'SCALARS scalars int '
            write(n_unit,'(A)') 'LOOKUP_TABLE default '
            DO II = 1, IIMX
                write(n_unit,'(I12)')NoB(II)
            END DO
        
            if (INP >= 1) then  !SHARP用。要素データを出力
                write(n_unit,'(A)') 'VECTORS vectors float'
                DO II = 1, IIMX
                    write(n_unit,*) UVW(II,1),UVW(II,2),UVW(II,3)
                END DO    
            end if
        close(n_unit)
        print*,'WRITEOUT:', FNAME
        
        !=======================================================================
    end subroutine writeout
        !**************************************************************************************


END MODULE grid_information
!==============================================================================================================
!境界面情報出力を追加(2021/04/08)
!飛沫計算用に特化させ、不要なサブルーチン削除(2021/04/27)
!**************************************************************************************************************
PROGRAM MAIN
!**************************************************************************************************************
    use grid_information
    IMPLICIT NONE
    character(8) :: d_start, d_stop
    character(10) :: t_start, t_stop
    character(99) :: FNAME, dir
    integer i
!==============================================================================================================
    call date_and_time(date = d_start, time = t_start)

    print*,'FILE NAME?'
    READ(5,'(A)') FNAME

    i = index(FNAME, "/", back=.true.)
    if(i <= 0) i = index(FNAME, "\", back=.true.)

    if(i > 0) dir = FNAME(1:i)

    INP = index(FNAME, '.inp') !INPファイルであれば自然数が返る。INPでないならゼロ。

    if(INP > 0) then
        call readINP(FNAME)
    else
        call readfile(FNAME)
    end if

    if(IIMX <= 0) then
        print*, 'ERROR_IIMX', IIMX
        stop
    end if

    call faceset
    call facecheck

    call boundaryset(dir)
    call nextcell(dir)   !必ずboundarysetの後にcall

    call writeout

    call date_and_time(date = d_stop, time = t_stop)
    print*, 'date = ', d_start, ' time = ', t_start
    print*, 'date = ', d_stop,  ' time = ', t_stop

END PROGRAM MAIN
