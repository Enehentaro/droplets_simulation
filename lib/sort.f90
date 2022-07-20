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
        integer switch

        contains

        procedure totalHeaplification!, partialHeaplification
        procedure get_featuredChildID, rebuild_tree
        
    end type

    public heap_sort
    public real2content

    contains

    !ヒープ木のコンストラクタ
    type(HeapTree) function HeapTree_(array)
        type(content_t), intent(in) :: array(:)

        HeapTree_%node = array

        call HeapTree_%totalHeaplification()

    end function

    !全体ヒープ化（親子の大小関係解決）メソッド
    subroutine totalHeaplification(self)
        class(HeapTree) self
        integer num_node
        integer parentID, child1ID, child2ID, featuredChildID

        num_node = size(self%node)
            
        do parentID = num_node/2, 1, -1 !子を持つノードに対してのみ下（葉の方）からループ

            child1ID = parentID*2
            child2ID = parentID*2 + 1

            featuredChildID = self%get_featuredChildID(child1ID, child2ID)

            if(self%node(parentID)%value > self%node(featuredChildID)%value) call swap_content(self%node, parentID, featuredChildID)

        end do

    end subroutine

    !入れ替えの起こった部分だけヒープ化（親子の大小関係解決）メソッド
    !こっちのほうが速いと思うが、現在バグってます
    ! subroutine partialHeaplification(self)
    !     class(HeapTree) self
    !     integer parentID, child1ID, child2ID, featuredChildID
    !     integer num_node

    !     num_node = size(self%node)

    !     parentID = 1

    !     do

    !         child1ID = parentID*2
    !         if(child1ID > num_node) exit    !第一子すら存在しなければループ終了
    !         child2ID = parentID*2 + 1

    !         featuredChildID = self%get_featuredChildID(child1ID, child2ID)

    !         if(self%node(parentID)%value > self%node(featuredChildID)%value) then
    !             call swap_content(self%node, parentID, featuredChildID)
    !             parentID = featuredChildID
    !         else
    !             exit    !入れ替えが起こらなければループ終了
    !         end if

    !     end do

    ! end subroutine

    !根ノードを除去し、ノード配列を左詰めにする。
    subroutine rebuild_tree(self)
        class(HeapTree) self

        self%node = self%node(2:)

    end subroutine

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

    subroutine heap_sort(array_origin, array_sorted)
        type(content_t), intent(in) :: array_origin(:)
        type(content_t), intent(out) :: array_sorted(:)
        type(HeapTree) heap_tree
        integer i
        integer arraySize

        arraySize = size(array_origin)
        if(size(array_sorted) /= arraySize) then
            print '("SORT ERROR")'
            error stop
        end if

        heap_tree = HeapTree_(array_origin)

        do i = 1, arraySize  !ソート後の配列に格納するループ
            array_sorted(i) = heap_tree%node(1)

            call heap_tree%rebuild_tree()

            call heap_tree%totalHeaplification()
            ! call heap_tree%partialHeaplification()

        end do

    end subroutine

    function real2content(real_array) result(content_array)
        real, intent(in) :: real_array(:)
        type(content_t), allocatable :: content_array(:)
        integer i

        allocate(content_array(size(real_array)))

        do i = 1, size(real_array)
            content_array(i)%originID = i
            content_array(i)%value = real_array(i)
        end do

    end function

end module sort_m