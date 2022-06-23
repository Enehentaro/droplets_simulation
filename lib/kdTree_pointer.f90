!井田の悪あがきです。これを実用化するつもりは今の所ありません。
module kdTree_pointer_m
    use sort_m
    implicit none
    private

    !ツリーノード構造体
    type node_in_kdTree
        type(node_in_kdTree), pointer :: parent=>null(), child1=>null(), child2=>null()
        integer :: cellID = 0
        type(content_t), allocatable :: array(:)
    end type

    !リスト内要素構造体
    !リスト内の要素としてツリーノードへのポインタを保持
    type content_in_list
        type(node_in_kdTree), pointer :: ptr=>null()
        type(content_in_list), pointer :: next=>null()
    end type

    !リストクラス
    type NodeList
        integer :: num_content = 0
        type(content_in_list), pointer :: head=>null(), tail=>null()
        contains
        procedure :: add => addToList
        procedure :: pop => popFromList
    end type

    !kdツリークラス
    type, public :: kdTree
        private
        type(node_in_kdTree), pointer :: root=>null(), end=>null()
        type(NodeList) leaf_list
        contains
        procedure :: addChildren => addChildrenToKdTree
        procedure :: print => print_KdTree
    end type

    public kdTree_

    contains

    !リストに要素を加えるメソッド
    subroutine addToList(self, node)
        class(NodeList) self
        type(node_in_kdTree), pointer :: node
        type(content_in_list), pointer :: content

        allocate(content)
        content%ptr => node

        if(self%num_content == 0) then
            self%head => content
            self%tail => content
        else
            self%tail%next => content
            self%tail => content
        end if

        self%num_content = self%num_content + 1

    end subroutine

    !リストの先頭から要素を取り出すメソッド
    function popFromList(self) result(node)
        class(NodeList) self
        type(node_in_kdTree), pointer :: node
        type(content_in_list), pointer :: new_head

        if(self%num_content >= 1) then
            node => self%head%ptr

            if(self%num_content == 1) then
                deallocate(self%head)
                self%head => null()
                self%tail => null()
            else
                new_head => self%head%next
                deallocate(self%head)
                self%head => new_head
            end if

            self%num_content = self%num_content - 1

        end if

    end function

    !kdツリー内の特定のノードを親として、子ノードを追加するメソッド
    subroutine addChildrenToKdTree(self, parent, child1, child2)
        class(kdTree) self
        type(node_in_kdTree), pointer :: parent, child1, child2

        if(associated(child1)) then
            parent%child1 => child1
            child1%parent => parent
            call self%leaf_list%add(child1) !新たなノードはリーフリストに追加
        end if

        if(associated(child2)) then
            parent%child2 => child2
            child2%parent => parent
            call self%leaf_list%add(child2) !新たなノードはリーフリストに追加
        end if

    end subroutine

    !kdツリーコンストラクタ
    type(kdTree) function kdTree_(array_origin)
        type(content_t), intent(in) :: array_origin(:)
        type(node_in_kdTree), pointer :: parent, child1, child2
        type(content_t), allocatable :: array(:), array_sorted(:)
        type(content_t), allocatable :: leftChildren(:), rightChildren(:)
        integer medianID

        allocate(kdTree_%root)
        call kdTree_%leaf_list%add(kdTree_%root)
        kdTree_%root%array = array_origin

        !リーフリストがなくなるまでループ
        do while(kdTree_%leaf_list%num_content >= 1)
            parent => kdTree_%leaf_list%pop()
            child1 => null()
            child2 => null()

            array = parent%array
            array_sorted = array !←配列サイズを揃えるため
            call heap_sort(array, array_sorted)
            medianID = size(array_sorted)/2 + 1
            parent%cellID = array_sorted(medianID)%originID

            leftChildren = array_sorted( : medianID-1 )
            rightChildren = array_sorted( medianID+1 : )

            if(size(leftChildren) >= 1) then
                allocate(child1)
                child1%array = leftChildren
            end if

            if(size(rightChildren) >= 1) then
                allocate(child2)
                child2%array = rightChildren
            end if

            call kdTree_%addChildren(parent, child1, child2)

        end do

    end function

    !kdツリーの内部を大まかに表示するメソッド
    subroutine print_KdTree(self)
        class(kdTree) self
        type(NodeList) open_list
        type(node_in_kdTree), pointer :: node
        integer i

        call open_list%add(self%root)

        do while(open_list%num_content >= 1)
            node => open_list%pop()

            print '("=======================================")'
            print*, node%cellID
            do i = 1, size(node%array)
                print '(i6, x, f20.16)', node%array(i)%originID, node%array(i)%value
            end do

            if(associated(node%child1)) then
                call open_list%add(node%child1)
            end if
    
            if(associated(node%child2)) then
                call open_list%add(node%child2)
            end if
            
        end do

    end subroutine

end module

program kdTree_test
    use kdTree_pointer_m
    use sort_m
    implicit none
    integer, parameter :: arraySize = 20
    type(content_t) array(arraySize)
    type(kdTree) kd_tree
    integer i

    call random_number(array(:)%value)
    do i = 1, arraySize
        array(i)%originID = i
    end do

    kd_tree = kdTree_(array)
    call kd_tree%print()

end program kdTree_test