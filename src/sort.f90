module sort_m
    implicit none
    private

    type, public :: content_t
        integer originID
        real value
        real coordinate(3)
    end type

    !ヒープ木クラス
    type HeapTree
        !実体は単なる配列だがツリー構造を表現している
        !要素 i に注目すると、親ノードは要素 i/2(小数切り捨て) であり、子ノードは要素 2i, 2i + 1 である
        type(content_t), allocatable :: node(:)
        integer switch

        contains

        procedure totalHeaplification, partialHeaplification
        procedure pop_from_root
        procedure get_featuredChildID
        
    end type

    public heap_sort

    contains

    !ヒープ木のコンストラクタ
    type(HeapTree) function HeapTree_(array, switch)
        type(content_t), intent(in) :: array(:)
        integer, intent(in) :: switch

        HeapTree_%node = array

        call HeapTree_%totalHeaplification()

        HeapTree_%switch = switch

    end function

    !全体ヒープ化（親子の大小関係解決）メソッド
    subroutine totalHeaplification(self)
        class(HeapTree) self
        integer num_node
        integer parentID, child1ID, child2ID, featuredChildID

        num_node = size(self%node)
        self%node(:)%value = self%node(:)%coordinate(self%switch)
            
        do parentID = num_node/2, 1, -1 !子を持つノードに対してのみ下（葉の方）からループ

            child1ID = parentID*2
            child2ID = parentID*2 + 1

            featuredChildID = self%get_featuredChildID(child1ID, child2ID)

            if(self%node(parentID)%value > self%node(featuredChildID)%value) call swap_content(self%node, parentID, featuredChildID)

        end do

    end subroutine

    !入れ替えの起こった部分だけヒープ化（親子の大小関係解決）メソッド
    subroutine partialHeaplification(self)
        class(HeapTree) self
        integer num_node, i
        type(content_t) endNode
        integer parentID, child1ID, child2ID, featuredChildID

        num_node = size(self%node)

        !末端ノードを根ノードに移動させ、全体をずらす
        endNode = self%node(num_node)
        do i = num_node, 2, -1
            self%node(i) = self%node(i-1)
        end do
        self%node(1) = endNode
        ! self%node = [self%node(num_node), self%node( : num_node-1 )]    !この書き方では遅い

        parentID = 1

        do

            child1ID = parentID*2
            if(child1ID > num_node) exit    !第一子すら存在しなければループ終了
            child2ID = parentID*2 + 1

            featuredChildID = self%get_featuredChildID(child1ID, child2ID)

            if(self%node(parentID)%value > self%node(featuredChildID)%value) then
                call swap_content(self%node, parentID, featuredChildID)
                parentID = featuredChildID
            else
                exit    !入れ替えが起こらなければループ終了
            end if

        end do

    end subroutine

    !ルートノードを返し、インスタンスからルートノードを除去する
    function pop_from_root(self) result(root)
        class(HeapTree) self
        type(content_t) root

        root = self%node(1)

        self%node = self%node(2:)

    end function

    subroutine swap_content(array, ID1, ID2)
        integer, intent(in) :: ID1, ID2
        type(content_t), intent(inout) :: array(:)
        type(content_t) temp !一時的に格納する変数

        temp = array(ID1)
        array(ID1) = array(ID2)
        array(ID2) = temp

    end subroutine

    function get_featuredChildID(self, child1ID, child2ID) result(featuredChildID)
        class(HeapTree) self
        integer, intent(in) :: child1ID, child2ID
        integer featuredChildID
        
        if(child2ID <= size(self%node)) then 
            !第2子が存在する（配列サイズ内）場合
            featuredChildID = get_smallerID(self%node(:)%value, child1ID, child2ID)

        else
            !第2子が存在しない（配列サイズ外）場合
            featuredChildID = child1ID  

        end if
    
    end function

    function get_smallerID(array, ID1, ID2) result(smaller_ID)
        real,intent(in) :: array(:)
        integer, intent(in) :: ID1, ID2
        integer smaller_ID

        if (array(ID1) < array(ID2)) then
            smaller_ID = ID1
        else
            smaller_ID = ID2
        end if
    
    end function

    ! switch = 1 なら x, switch = 2 なら y, switch = 3 なら z
    subroutine heap_sort(array_origin, switch, array_sorted)
        type(content_t), intent(in) :: array_origin(:)
        integer, intent(in) :: switch
        type(content_t), intent(out), allocatable :: array_sorted(:)
        type(HeapTree) heap_tree
        integer i
        integer arraySize

        arraySize = size(array_origin)
        allocate(array_sorted(arraySize))
        ! if(size(array_sorted) /= arraySize) then
        !     print '("SORT ERROR")'
        !     stop
        ! end if

        heap_tree = HeapTree_(array_origin, switch)

        do i = 1, arraySize  !ソート後の配列に格納するループ
            array_sorted(i) = heap_tree%pop_from_root()
            ! call heap_tree%totalHeaplification()
            call heap_tree%partialHeaplification()
        end do

    end subroutine

end module sort_m