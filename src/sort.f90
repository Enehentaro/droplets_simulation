module sort_m
    implicit none
    private

    type, public :: content_t
        integer originID
        real value
    end type

    !ヒープ木クラス
    type HeapTree
        !実体は単なる配列だがツリー構造を表現している
        !要素 i に注目すると、親ノードは要素 i/2(小数切り捨て) であり、子ノードは要素 2i, 2i + 1 である
        type(content_t), allocatable :: node(:)

        contains

        procedure heaplification, pop_from_root
        
    end type

    public heap_sort

    contains

    !ヒープ木のコンストラクタ
    type(HeapTree) function HeapTree_(array)
        type(content_t), intent(in) :: array(:)

        HeapTree_%node = array
        call HeapTree_%heaplification()

    end function

    !ヒープ化（親子の大小関係解決）メソッド
    subroutine heaplification(self)
        class(HeapTree) self
        integer num_node
        integer parentID, child1ID, child2ID, featuredChildID

        num_node = size(self%node)
            
        do parentID = num_node/2, 1, -1 !子を持つノードに対してのみ下からループ

            child1ID = parentID*2
            child2ID = parentID*2 + 1

            if(child2ID <= num_node) then !第2子が存在する（配列サイズ内）の場合

                featuredChildID = get_smallerID(self%node(:)%value, child1ID, child2ID)

            else

                featuredChildID = child1ID  !第2子が存在しない（配列サイズ外）の場合

            end if

            if(self%node(parentID)%value > self%node(featuredChildID)%value) call swap_content(self%node, parentID, featuredChildID)

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
    
    subroutine heap_sort(array_origin, array_sorted)
        type(content_t), intent(in) :: array_origin(:)
        type(content_t), intent(out) :: array_sorted(:)
        type(HeapTree) heap_tree
        integer i
        integer arraySize

        arraySize = size(array_origin)
        if(size(array_sorted) /= arraySize) then
            print '("SORT ERROR")'
            stop
        end if

        heap_tree = HeapTree_(array_origin)

        do i = 1, arraySize  !ソート後の配列に格納するループ
            array_sorted(i) = heap_tree%pop_from_root()
            call heap_tree%heaplification()
        end do

    end subroutine

end module sort_m