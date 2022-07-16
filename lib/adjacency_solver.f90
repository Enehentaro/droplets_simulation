MODULE adjacencySolver_m
    ! use unstructuredGrid_m
    IMPLICIT NONE
    private
    ! integer num_nodes, num_cells, num_tetras, num_prisms, num_pyramids
    ! INTEGER num_halfFace

    type halfFace_t
        integer :: nodeID(4) = 0, pairID = 0, ID_sum = 0, ownerID(2) = -1
    end type halfFace_t

    type AdjacencySolver
        type(halfFace_t), allocatable :: halfFaceArray(:)
        integer num_BoundFace
        contains
        procedure set_halfFaceArray,  check_halfFace, find_boundFaceInformation, find_adjacentCellID
    end type

    public solve_BoundaryAndAdjacency

    contains

    subroutine set_halfFaceArray(self, cellVertices)
        class(AdjacencySolver) self
        integer, intent(in) :: cellVertices(:,:)
        INTEGER II,JJJ, j, n, JJJMX, num_halfFace, num_cell
        integer, allocatable :: n_type(:)
        integer, parameter :: num_halfFace_perCELL(3) = [4,5,5]
        integer, parameter :: IDtrans(4,5,3) = reshape([ &
                                1,2,3,0, 2,3,4,0, 3,4,1,0, 4,1,2,0, 0,0,0,0,&                       !テトラ
                                1,2,3,0, 4,5,6,0, 1,2,4,5, 2,3,5,6, 3,1,6,4,&                       !プリズム
                                5,1,2,0, 5,2,3,0, 5,3,4,0, 5,4,1,0, 1,2,3,4 ], shape(IDtrans))      !ピラミッド
                                                        
        num_cell = size(cellVertices, dim=2)    !cellVerticesには、(6ｘセル数)の配列が入ってきている想定
        if(num_cell <= 0) then
            print*, 'ERROR_num_cells', num_cell
            error stop
        end if

        allocate(n_type(num_cell))
        DO II = 1, num_cell
            select case(count(cellVertices(:,II)==0))
                case(2) !ゼロの数が2個：頂点数が4個：テトラ
                    n_type(II) = 1
                case(0) !ゼロの数が0個：頂点数が6個：プリズム
                    n_type(II) = 2
                case(1) !ゼロの数が1個：頂点数が5個：ピラミッド
                    n_type(II) = 3
                case default
                    print*, '**CEll VERTICES ERROR**'
                    error stop
            end select
        END DO

        num_halfFace = count(n_type==1)*4 + count(n_type==2)*5 + count(n_type==3)*5   !半面数：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
    
        allocate(self%halfFaceArray(num_halfFace))
            
        JJJ = 0
        
        DO II = 1, num_cell


            do j = 1, num_halfFace_perCELL(n_type(II))
                JJJ = JJJ + 1
                do n = 1, 3
                    ! print*, JJJ, n, j, n_type, II
                    self%halfFaceArray(JJJ)%nodeID(n) = cellVertices(IDtrans(n, j, n_type(II)), II)

                end do
                if(IDtrans(4, j, n_type(II)) > 0) then
                    self%halfFaceArray(JJJ)%nodeID(4) = cellVertices(IDtrans(4, j, n_type(II)), II) !四角形面なら4点目を代入

                end if
                self%halfFaceArray(JJJ)%ownerID(1) = II
                self%halfFaceArray(JJJ)%ID_sum = sum(self%halfFaceArray(JJJ)%nodeID(:))

            end do
            
        END DO    
              
        JJJMX = JJJ
        
        if (JJJMX == num_halfFace) then
            print*, 'JJJMX / num_halfFace =', JJJMX, '/', num_halfFace
        else
            print*, 'JJJMX_ERROR:', JJJMX, num_halfFace
            error stop
        end if
            
    end subroutine
    
    subroutine check_halfFace(self)
        use terminalControler_m
        class(AdjacencySolver) self
        INTEGER match, width, numNode, faceID, groupID,num_group, maxID_sum, num_halfFace, num_BoundFaces
        integer checkCounter, faceCounter, i,j, faceID1,faceID2, num_face, k,l
        integer, allocatable :: faceID_array(:)
        type faceGroup_t
            integer, allocatable :: faceID(:)
        end type faceGroup_t
        type(faceGroup_t), allocatable :: faceGroup(:)
        ! real time1, time2

        ! call cpu_time(time1)

        print*,'START-FACE CHECK!' !同一面の探索
        
        num_halfFace = size(self%halfFaceArray)
        num_group = num_halfFace/5000 + 1   !面グループ数（1グループ数に約5000枚面が入るようにする）（この値は経験則）
    
        allocate(faceGroup(num_group))
        allocate(faceID_array(num_halfFace), source=0)

        maxID_sum = maxval(self%halfFaceArray(:)%ID_sum)
        width = maxID_sum/num_group + 1   !1グループの幅（最大節点番号和を面グループ数で割る）
          
        checkCounter = 0
        call set_formatTC('("DIVIDE halfFace [ #group : ",i6," / ",i6," ]")')
        do groupID = 1, num_group  !面をグループに分ける
            call print_progress([groupID, num_group])
            faceID_array(:) = 0
            faceCounter = 1
            do faceID = 1, num_halfFace
                if((self%halfFaceArray(faceID)%ID_sum > (groupID-1)*width) &
                    .and.(self%halfFaceArray(faceID)%ID_sum <= groupID*width)) then

                    faceID_array(faceCounter) = faceID
                    faceCounter = faceCounter + 1
                    checkCounter = checkCounter + 1
                    
                end if
            end do
            faceGroup(groupID)%faceID = faceID_array(1 : faceCounter-1)
        end do

        if(checkCounter /= num_halfFace) then
            print*, 'faceCounter_ERROR:', checkCounter, num_halfFace
            error stop
        end if
          
        num_BoundFaces = 0
        call set_formatTC('("CHECK halfFace [ #group : ",i6," / ",i6," ]")')
        !$omp parallel do private(faceID1,faceID2, match, k,l, num_face,numNode) reduction(+:num_BoundFaces)
        do groupID = 1, num_group
            call print_progress([groupID, num_group])
            num_face = size(faceGroup(groupID)%faceID)

            face1 : do i = 1, num_face
                faceID1 = faceGroup(groupID)%faceID(i)
                if(self%halfFaceArray(faceID1)%ownerID(2) >= 0) cycle face1 !面共有セル探索が済んでいる場合スキップ
                if(self%halfFaceArray(faceID1)%nodeID(4) <= 0) then
                    numNode = 3   !三角形面
                else
                    numNode = 4   !四角形面
                end if
        
                face2 :do j = i+1, num_face
                    faceID2 = faceGroup(groupID)%faceID(j)
                    if(self%halfFaceArray(faceID2)%ownerID(2) >= 0) cycle face2 !面共有セル探索が済んでいる場合スキップ
                    if((numNode==3).and.(self%halfFaceArray(faceID2)%nodeID(4) > 0)) cycle face2 !三角形面を注目中に四角形面が現れればスキップ
            
                    match = 0
                    do k = 1, numNode
                        do l = 1, numNode
                            if(self%halfFaceArray(faceID1)%nodeID(k) == self%halfFaceArray(faceID2)%nodeID(l)) match = match + 1  !点IDが一致すればカウント
                        end do
                    end do
            
                    if(match < numNode) cycle face2 !節点数と同じ回数一致しなければスキップ（必要十分条件）
            
                    !ここまでくれば同一の2面発見
            
                    self%halfFaceArray(faceID1)%ownerID(2) = self%halfFaceArray(faceID2)%ownerID(1)
                    self%halfFaceArray(faceID2)%ownerID(2) = self%halfFaceArray(faceID1)%ownerID(1)
            
                    cycle face1 !共有セルが見つかったので次の面へ
        
                end do face2
        
                self%halfFaceArray(faceID1)%ownerID(2) = 0 !共有セルが見つからなかった:境界面
                num_BoundFaces = num_BoundFaces + 1  !境界面カウント
        
            end do face1
    
        end do
        !$omp end parallel do
        
        print*,'# Boundary Face =', num_BoundFaces

        self%num_BoundFace = num_BoundFaces

        print*,'END-FACE CHECK!'

        ! call cpu_time(time2)
        ! print*, time2 - time1
        ! error stop

    END subroutine check_halfFace

    subroutine find_boundFaceInformation(self, cellBoundFaces, boundFaceVertices)
        class(AdjacencySolver) self
        INTEGER II,JJJ,JB
        integer cellBoundFaces(:,:)
        integer, allocatable :: NoB(:)
        integer, allocatable, intent(out) :: boundFaceVertices(:,:)

        allocate(NoB(size(cellBoundFaces, dim=2)), source=0)
        allocate(boundFaceVertices(3, self%num_BoundFace), source=0)

        JB = 0
        DO JJJ = 1, size(self%halfFaceArray)
            IF (self%halfFaceArray(JJJ)%ownerID(2) /= 0) cycle   !境界面以外はスルー
            JB = JB +1
            II = self%halfFaceArray(JJJ)%ownerID(1) !JJが属する要素番号
            NoB(II) = NoB(II) +1 !IIが所有する境界面数カウント
            cellBoundFaces(NoB(II),II) = JB !IIが所有する境界面番号

            boundFaceVertices(1:3, JB) = self%halfFaceArray(JJJ)%nodeID(1:3)

        END DO
       
    end subroutine
        
    subroutine find_adjacentCellID(self, adjacentCellArray)
        class(AdjacencySolver) self
        integer II, JJJ, adjacentCellArray(:,:)
        integer, allocatable :: num_adjacent(:)

        allocate(num_adjacent(size(adjacentCellArray)), source=0)

        do JJJ = 1, size(self%halfFaceArray)
            if (self%halfFaceArray(JJJ)%ownerID(2) <= 0) cycle !境界面はスルー

            II = self%halfFaceArray(JJJ)%ownerID(1)
            num_adjacent(II) = num_adjacent(II) + 1
            adjacentCellArray(num_adjacent(II), II) = self%halfFaceArray(JJJ)%ownerID(2)

        end do

    end subroutine

    subroutine solve_BoundaryAndAdjacency(cellVertices, cellBoundFaces, boundFaceVertices, adjacentCellArray)
        integer, intent(in) :: cellVertices(:,:)
        integer cellBoundFaces(:,:), adjacentCellArray(:,:)
        integer, allocatable, intent(out) :: boundFaceVertices(:,:)
        type(AdjacencySolver) AS

        print*, '~ SolvingAdjacency is Required. ~'

        call AS%set_halfFaceArray(cellVertices)    !面情報のセッティング

        call AS%check_halfFace()  !同一面のチェック

        call AS%find_boundFaceInformation(cellBoundFaces, boundFaceVertices)   !境界面情報

        call AS%find_adjacentCellID(adjacentCellArray)   !セル隣接情報

    end subroutine

END MODULE adjacencySolver_m