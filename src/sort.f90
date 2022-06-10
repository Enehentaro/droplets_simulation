module sort_m
    implicit none
    ! private

    type content
        integer originID
        real axis
    end type content

    contains

    subroutine swap_content(array, ID1, ID2)
        integer, intent(in) :: ID1, ID2
        type(content), intent(inout) :: array(:)
        type(content) temp !一時的に格納する変数

        temp = array(ID1)
        array(ID1) = array(ID2)
        array(ID2) = temp

    end subroutine

    function get_largerID(array, ID1, ID2) result(larger_ID)
        real,intent(in) :: array(:)
        integer, intent(in) :: ID1, ID2
        integer larger_ID

        if (array(ID1) > array(ID2)) then
            larger_ID = ID1
        else
            larger_ID = ID2
        end if
    
    end function
    
    subroutine heap_sort(array_origin, ID_origin, array_sorted)
        real, intent(in) :: array_origin(:)
        type(content), intent(out), allocatable :: array_sorted(:)
        integer, intent(in) :: ID_origin(:)
        type(content), allocatable :: pre_array(:)

        type(content), allocatable :: temp_array(:)  !配列の要素を減らすための一時的な配列
        type(content) parent, child, larger_child, child_1st, child_2nd
        integer i, j, node, num_node, num_node_heap
        integer parentID, largerChildID

        num_node = size(array_origin)
        num_node_heap = num_node    !初めは要素全体のノード数とヒープ構造のノード数が同じ
        allocate(pre_array(num_node))
        do i = 1, num_node
            pre_array(i)%originID = ID_origin(i)
            pre_array(i)%axis = array_origin(i)
        end do
        allocate(array_sorted(num_node))

        do i = 1, num_node  !ソート後の配列に格納するループ

            if (mod(num_node_heap,2) == 0) then  !要素数の偶奇判定(偶数のときだけ特別な処理)
                parent = pre_array(num_node_heap/2)
                child = pre_array(num_node_heap)

                if(parent%axis < child%axis) then !要素数が偶数のとき末端のノードだけ2分木にならないのでその処理
                    call swap_content(pre_array, num_node_heap/2, num_node_heap)
                end if

                num_node_heap = num_node_heap-1   !末端の1分木の分だけ比較に使用するヒープ構造の要素数を減らす
            
            end if
                
            do node = num_node_heap, 3, -2   !ヒープソート一回分のループ

                parentID = int(node/2)

                largerChildID = get_largerID(pre_array(:)%axis, node-1, node)

                if(pre_array(parentID)%axis < pre_array(largerChildID)%axis) then
                    call swap_content(pre_array, parentID, largerChildID)
                end if

            end do

            array_sorted(i) = pre_array(1)
            allocate(temp_array(size(pre_array)-1))

            do j = 2, size(pre_array)  !入力配列の最上親ノードを抜き取って1つずらした配列作成
                temp_array(j-1) = pre_array(j)
            end do

            deallocate(pre_array)
            allocate(pre_array(num_node-i))

            pre_array = temp_array

            deallocate(temp_array)

            num_node_heap = size(pre_array)

        end do

    end subroutine

end module sort_m