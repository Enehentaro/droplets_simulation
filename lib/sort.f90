module sort_m
    implicit none
    private

    !>コンテンツ構造体
    !>実数とIDをメンバに持つ
    type, public :: content_t
        integer originID
        real value
    end type

    !>ヒープ木クラス
    !>実体は単なる配列だがツリー構造を表現している
    !>要素 i に注目すると、親ノードは要素 i/2(小数切り捨て) であり、子ノードは要素 2i, 2i + 1 である
    type HeapTree
        type(content_t), allocatable :: node(:)

        contains

        procedure totalHeaplification!, partialHeaplification
        procedure rebuild_tree!, get_featuredChildID
        
    end type

    public heap_sort
    public real2content
    public merge_sort

    contains

    !>ヒープ木のコンストラクタ
    type(HeapTree) function HeapTree_(array)
        type(content_t), intent(in) :: array(:)

        HeapTree_%node = array

        call HeapTree_%totalHeaplification()

    end function

    !>全体ヒープ化（親子の大小関係解決）メソッド
    subroutine totalHeaplification(self)
        class(HeapTree) self
        integer num_node
        integer parentID, child1ID, child2ID!, featuredChildID

        num_node = size(self%node)
            
        do parentID = num_node/2, 1, -1 !子を持つノードに対してのみ下（葉の方）からループ

            child1ID = parentID*2
            child2ID = parentID*2 + 1

            ! featuredChildID = self%get_featuredChildID(child1ID, child2ID)
            !これ使わんほうが速い

            if(self%node(parentID)%value > self%node(child1ID)%value) call swap_content(self%node, parentID, child1ID)

            if(child2ID <= num_node) then
                if(self%node(parentID)%value > self%node(child2ID)%value) call swap_content(self%node, parentID, child2ID)
            end if

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

    subroutine rebuild_tree(self)
        !!根ノードを除去し、ノード配列を左詰めにする。
        !!このとき、配列のサイズが一つ減る。
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

    ! function get_featuredChildID(self, child1ID, child2ID) result(featuredChildID)
    !     class(HeapTree) self
    !     integer, intent(in) :: child1ID, child2ID
    !     integer featuredChildID
        
    !     if(child2ID <= size(self%node)) then 
    !         !第2子が存在する（配列サイズ内）場合
    !         !値の小さい方を返す
    !         if (self%node(child1ID)%value < self%node(child2ID)%value) then
    !             featuredChildID = child1ID
    !         else
    !             featuredChildID = child2ID
    !         end if

    !         !self%node(:)%valueの臨時配列生成に時間がかかってたぽい
    !         ! featuredChildID = get_smallerID(self%node(:)%value, child1ID, child2ID)

    !     else
    !         !第2子が存在しない（配列サイズ外）場合
    !         featuredChildID = child1ID  

    !     end if
    
    ! end function

    ! function get_smallerID(array, ID1, ID2) result(smaller_ID)
    !     real,intent(in) :: array(:)
    !     integer, intent(in) :: ID1, ID2
    !     integer smaller_ID

    !     if (array(ID1) < array(ID2)) then
    !         smaller_ID = ID1
    !     else
    !         smaller_ID = ID2
    !     end if
    
    ! end function

    !>ヒープソート
    function heap_sort(array_origin) result(array_sorted)
        type(content_t), intent(in) :: array_origin(:)
        type(content_t) array_sorted(size(array_origin))
        type(HeapTree) heap_tree
        integer i

        heap_tree = HeapTree_(array_origin)

        do i = 1, size(array_origin)    !ソート後の配列に格納するループ
            array_sorted(i) = heap_tree%node(1)

            call heap_tree%rebuild_tree()

            call heap_tree%totalHeaplification()
            ! call heap_tree%partialHeaplification()

        end do

    end function

    !>実数型配列をコンテンツ配列に変換する
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

    recursive subroutine merge_sort(array, left_end, right_end)
        type(content_t), intent(inout) :: array(:)
        integer, intent(in) :: left_end, right_end
        type(content_t), allocatable :: work_array(:)
        integer mid, i,j,k

        allocate(work_array(size(array)))
        
        if(left_end < right_end) then
            mid = int((left_end + right_end)/2)
            call merge_sort(array, left_end, mid)
            call merge_sort(array, mid+1, right_end)
            
            do i = mid, left_end, -1
                work_array(i) = array(i)           
            end do

            do j = mid+1, right_end
                work_array(right_end-(j-(mid+1))) = array(j)
            end do

            i = left_end
            j = right_end

            do k = left_end, right_end
                if(work_array(i)%value < work_array(j)%value) then
                    array(k) = work_array(i)
                    i = i + 1
                else
                    array(k) = work_array(j)
                    j = j - 1
                end if
            end do

        end if

    end subroutine

end module sort_m