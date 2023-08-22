MODULE adjacencySolver_m
    !!セルの隣接関係解決モジュール
    IMPLICIT NONE
    private

    integer, parameter, public :: None = -1
        !!整数配列中の欠損値の表現

    integer, parameter :: max_adjacent=5  !想定される最大隣接セル数（プリズムは最大５つのセルと隣接）
    integer, parameter :: max_boundFace=4  !想定される最大境界面数（プリズムの５つの面のうち、すべて境界面はありえないので４）

    !>ハーフフェイス構造体
    type halfFace_t
        integer :: nodeID(4) = None, pairID = None, ID_sum = 0, ownerCellID(2) = None
    end type

    public solve_BoundaryAndAdjacency

    contains

    function get_halfFaceArray(cellVertices) result(halfFaceArray)
        !!各セルごとにハーフフェイスをカウントし、配列に格納
        integer, intent(in) :: cellVertices(:,:)
        type(halfFace_t), allocatable :: halfFaceArray(:)
        INTEGER II,JJJ, j, n, JJJMX, num_halfFace, num_cell, vertexID
        character(5), allocatable :: cellType(:)
        character(5), parameter :: tetra='tetra', prism='prism', pyramid='pyrmd'    !キーワード
        integer, allocatable :: halfFaceVerticesID(:,:)

        !以下、各ハーフフェイス頂点の順番を表す配列。
        !例えば、テトラは４つの頂点から構成されており、ハーフフェイスも４つ。
        !１つ目のハーフフェイスは、４つの頂点のうち 1,2,3 番目のものを頂点とする。
        !以降、 2,3,4, 3,4,1, 4,1,2 番目を頂点とし、４つのハーフフェイスの定義が完了する。
        integer, parameter :: halfFaceVerticesID_inTetra(3,4) = reshape(&
            [1,2,3, 2,3,4, 3,4,1, 4,1,2], shape(halfFaceVerticesID_inTetra))

        integer, parameter :: halfFaceVerticesID_inPrism(4,5) = reshape(&
            [1,2,3,None, 4,5,6,None, 1,2,4,5, 2,3,5,6, 3,1,6,4], shape(halfFaceVerticesID_inPrism))

        integer, parameter :: halfFaceVerticesID_inPyramid(4,5) = reshape(&
            [5,1,2,None, 5,2,3,None, 5,3,4,None, 5,4,1,None, 1,2,3,4], shape(halfFaceVerticesID_inPyramid))
                           
            
        num_cell = size(cellVertices, dim=2)    !cellVerticesには、(6ｘセル数)の配列が入ってきている想定
        if(num_cell <= 0) then
            print*, 'ERROR_num_cells', num_cell
            error stop
        end if

        allocate(cellType(num_cell))
        DO II = 1, num_cell
            select case(count(cellVertices(:,II)/=None))    !頂点数（配列内のNoneでない要素数）で場合分け
            case(4) !頂点数が4個：テトラ
                cellType(II) = tetra

            case(6) !頂点数が6個：プリズム
                cellType(II) = prism

            case(5) !頂点数が5個：ピラミッド
                cellType(II) = pyramid

            case default
                print*, '**CEll VERTICES ERROR**'
                error stop
            end select
        END DO

        num_halfFace = count(cellType==tetra)*4 + count(cellType==prism)*5 + count(cellType==pyramid)*5
            !半面数：テトラ数×4 (+プリズム数×5 +ピラミッド数×5)
    
        allocate(halfFaceArray(num_halfFace))
            
        JJJ = 0
        
        DO II = 1, num_cell
            select case(cellType(II))
            case(tetra)
                halfFaceVerticesID = halfFaceVerticesID_inTetra
            case(prism)
                halfFaceVerticesID = halfFaceVerticesID_inPrism
            case(pyramid)
                halfFaceVerticesID = halfFaceVerticesID_inPyramid
            end select

            !halfFaceVerticesIDは２次元配列で、（頂点、ハーフフェイス）

            do j = 1, size(halfFaceVerticesID, dim=2)!ハーフフェイスの数だけループ
                JJJ = JJJ + 1   !ハーフフェイスの絶対IDをカウント

                do n = 1, size(halfFaceVerticesID, dim=1)!あるハーフフェイスにおける頂点数だけループ
                    vertexID = halfFaceVerticesID(n,j)
                    if(vertexID/=None) halfFaceArray(JJJ)%nodeID(n) = cellVertices(vertexID, II)
                end do

                halfFaceArray(JJJ)%ownerCellID(1) = II  !現在はセルIIに着目しているので、自明ながらオーナーセルもII

                halfFaceArray(JJJ)%ID_sum = sum(halfFaceArray(JJJ)%nodeID(:))   !IDの和を格納（あとで使う）

            end do
            
        END DO    
              
        JJJMX = JJJ
        
        if (JJJMX == num_halfFace) then
            print*, 'JJJMX / num_halfFace =', JJJMX, '/', num_halfFace
        else
            print*, 'JJJMX_ERROR:', JJJMX, num_halfFace
            error stop
        end if
            
    end function
    
    subroutine check_halfFace(halfFaceArray, num_BoundFaces)
        !!各ハーフフェイスに対して相方を探す
        !!相方がみつかれば、相方のセルと隣接していることがわかる
        !!相方のいないハーフフェイスは境界面
        use terminalControler_m
        type(halfFace_t), intent(inout) :: halfFaceArray(:)
        integer, intent(out) :: num_BoundFaces
        INTEGER match, width, numNode, faceID, groupID,num_group, maxID_sum, num_halfFace
        integer checkCounter, faceCounter, i,j, faceID1,faceID2, num_face, k,l
        integer, allocatable :: faceID_array(:)
        type faceGroup_t
            integer, allocatable :: faceID(:)
        end type faceGroup_t
        type(faceGroup_t), allocatable :: faceGroup(:)
        ! real time1, time2

        ! call cpu_time(time1)

        print*,'START-FACE CHECK!' !同一面の探索
        
        num_halfFace = size(halfFaceArray)
        num_group = num_halfFace/5000 + 1   !面グループ数（1グループ数に約5000枚面が入るようにする）（この値は経験則）
    
        allocate(faceGroup(num_group))
        allocate(faceID_array(num_halfFace), source=0)

        maxID_sum = maxval(halfFaceArray(:)%ID_sum)
        width = maxID_sum/num_group + 1   !1グループの幅（最大節点番号和を面グループ数で割る）
          
        checkCounter = 0
        call set_formatTC('("DIVIDE halfFace [ #group : ",i6," / ",i6," ]")')
        do groupID = 1, num_group  !面をグループに分ける
            call print_progress([groupID, num_group])
            faceID_array(:) = 0
            faceCounter = 1
            do faceID = 1, num_halfFace
                if((halfFaceArray(faceID)%ID_sum > (groupID-1)*width) &
                    .and.(halfFaceArray(faceID)%ID_sum <= groupID*width)) then

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
        !!$omp parallel do private(faceID1,faceID2, match, k,l, num_face,numNode) reduction(+:num_BoundFaces)
        do groupID = 1, num_group
            call print_progress([groupID, num_group])
            num_face = size(faceGroup(groupID)%faceID)

            face1 : do i = 1, num_face
                faceID1 = faceGroup(groupID)%faceID(i)
                if(halfFaceArray(faceID1)%ownerCellID(2) /= None) cycle face1 !面共有セル探索が済んでいる場合スキップ
                if(halfFaceArray(faceID1)%nodeID(4) == None) then
                    numNode = 3   !三角形面
                else
                    numNode = 4   !四角形面
                end if
        
                face2 :do j = i+1, num_face
                    faceID2 = faceGroup(groupID)%faceID(j)
                    if(halfFaceArray(faceID2)%ownerCellID(2) /= None) cycle face2 !面共有セル探索が済んでいる場合スキップ
                    if((numNode == 3) .and. (halfFaceArray(faceID2)%nodeID(4) /= None)) cycle face2 !三角形面を注目中に四角形面が現れればスキップ
            
                    match = 0
                    do k = 1, numNode
                        do l = 1, numNode
                            if(halfFaceArray(faceID1)%nodeID(k) == halfFaceArray(faceID2)%nodeID(l)) match = match + 1  !点IDが一致すればカウント
                        end do
                    end do
            
                    if(match < numNode) cycle face2 !節点数と同じ回数一致しなければスキップ（必要十分条件）
            
                    !ここまでくれば同一の2面発見
            
                    halfFaceArray(faceID1)%ownerCellID(2) = halfFaceArray(faceID2)%ownerCellID(1)
                    halfFaceArray(faceID2)%ownerCellID(2) = halfFaceArray(faceID1)%ownerCellID(1)
            
                    cycle face1 !共有セルが見つかったので次の面へ
        
                end do face2
        
                ! halfFaceArray(faceID1)%ownerCellID(2) = 0 !共有セルが見つからなかった:境界面
                num_BoundFaces = num_BoundFaces + 1  !境界面カウント
        
            end do face1
    
        end do
        !!$omp end parallel do
        
        print*,'# Boundary Face =', num_BoundFaces

        print*,'END-FACE CHECK!'

        ! call cpu_time(time2)
        ! print*, time2 - time1
        ! error stop

    END subroutine check_halfFace

    subroutine find_boundFaceInformation(halfFaceArray, num_boundFace, cellBoundFaces, triangleBoundFaceVertices)
        !!境界面のみを取り出し、配列に格納
        !!四角形面であったとしても、最初の3点だけ抽出し、三角形面にして返す
        type(halfFace_t), intent(in) :: halfFaceArray(:)
        integer, intent(in) :: num_boundFace
        INTEGER II,JJJ,JB
        integer, intent(out) :: cellBoundFaces(:,:)
        integer, allocatable :: NoB(:)
        integer, allocatable, intent(out) :: triangleBoundFaceVertices(:,:)

        allocate(NoB(size(cellBoundFaces, dim=2)), source=0)
        allocate(triangleBoundFaceVertices(3, num_BoundFace))

        JB = 0
        DO JJJ = 1, size(halfFaceArray)
            if (halfFaceArray(JJJ)%ownerCellID(2) /= None) cycle   !第二オーナーが存在：境界面ではない：スルー

            JB = JB +1
            II = halfFaceArray(JJJ)%ownerCellID(1) !JJが属する要素番号
            NoB(II) = NoB(II) +1 !IIが所有する境界面数カウント
            cellBoundFaces(NoB(II),II) = JB !IIが所有する境界面番号

            triangleBoundFaceVertices(1:3, JB) = halfFaceArray(JJJ)%nodeID(1:3) !四角形面であったとしても、最初の3点だけ抽出

        END DO
       
    end subroutine
        
    subroutine find_adjacentCellID(halfFaceArray, adjacentCellIDArray)
        !!隣接関係を配列に格納
        type(halfFace_t), intent(in) :: halfFaceArray(:)
        integer, intent(out) :: adjacentCellIDArray(:,:)
        integer II, JJJ
        integer, allocatable :: num_adjacent(:)

        allocate(num_adjacent(size(adjacentCellIDArray)), source=0)

        do JJJ = 1, size(halfFaceArray)
            if (halfFaceArray(JJJ)%ownerCellID(2) == None) cycle !第２オーナーが未発見のハーフフェイス：境界面：スルー

            II = halfFaceArray(JJJ)%ownerCellID(1)
            num_adjacent(II) = num_adjacent(II) + 1
            adjacentCellIDArray(num_adjacent(II), II) = halfFaceArray(JJJ)%ownerCellID(2)

        end do

    end subroutine

    subroutine solve_BoundaryAndAdjacency(cellVertices, cellBoundFaces, triangleBoundFaceVertices, adjacentCellIDArray)
        !!境界面と隣接関係を、それぞれ配列に格納

        integer, intent(in) :: cellVertices(:,:)
            !!セルの頂点ID配列（頂点ID,セルID）

        integer, allocatable, intent(out) :: cellBoundFaces(:,:)
            !!セルの境界面ID配列（境界面ID,セルID）

        integer, allocatable, intent(out) :: triangleBoundFaceVertices(:,:)
            !!境界面の頂点ID配列（頂点ID,境界面ID）

        integer, allocatable, intent(out) :: adjacentCellIDArray(:,:)
            !!セルの隣接セルID配列（隣接セルID,セルID）

        type(halfFace_t), allocatable :: halfFaceArray(:)
        integer num_cell, num_BoundFace

        print*, '~ SolvingAdjacency is Required. ~'

        halfFaceArray = get_halfFaceArray(cellVertices)    !ハーフフェイスのセッティング

        call check_halfFace(halfFaceArray, num_BoundFace)  !同一面のチェック

        num_cell = size(cellVertices, dim=2)
        allocate(cellBoundFaces(max_boundFace, num_cell), source=None)
        allocate(adjacentCellIDArray(max_adjacent, num_cell), source=None)
        call find_boundFaceInformation(halfFaceArray, num_BoundFace, cellBoundFaces, triangleBoundFaceVertices)   !境界面情報

        call find_adjacentCellID(halfFaceArray, adjacentCellIDArray)   !セル隣接情報

    end subroutine

END MODULE adjacencySolver_m