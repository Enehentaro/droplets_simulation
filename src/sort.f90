module sort_m
    implicit none
    private

    contains

    subroutine swap_node(parent, child)
        real, intent(inout) :: parent, child
        real temp   !一時的に格納する変数

        temp = parent
        parent = child
        child = temp

    end subroutine

    function compare_children(child_1st, child_2nd) result(child)
        real, intent(in) :: child_1st, child_2nd
        real child

        if (child_1st < child_2nd) then
            child = child_2nd
        else
            child = child_1st
        end if
    
    end function
    
    subroutine heap_sort(array, array_sorted)
        real, intent(inout), allocatable :: array(:)
        real, intent(out), allocatable :: array_sorted(:)
        real, allocatable :: temp_array(:)  !配列の要素を減らすための一時的な配列
        real parent, child, child_1st, child_2nd
        integer i, j, node, num_node, num_node_heap

        num_node = size(array)
        num_node_heap = num_node    !初めは要素全体のノード数とヒープ構造のノード数が同じ
        allocate(array_sorted(num_node))
        do i = 1, num_node  !ソート後の配列に格納するループ

            if (mod(num_node_heap,2) == 0) then  !要素数の偶奇判定(偶数のときだけ特別な処理)
                parent = array(num_node_heap/2)
                child = array(num_node_heap)

                if(parent < child) then !要素数が偶数のとき末端のノードだけ2分木にならないのでその処理
                    call swap_node(parent, child)
                    array(num_node_heap/2) = parent
                    array(num_node_heap) = child
                end if

                num_node_heap = num_node_heap-1   !末端の1分木の分だけ比較に使用するヒープ構造の要素数を減らす
            
            end if
                
            do node = num_node_heap, 3, -2   !ヒープソート一回分のループ
                parent = array(int(node/2))
                child_1st = array(node-1)
                child_2nd = array(node)
                child = compare_children(child_1st, child_2nd)
                if(parent < child) then
                    call swap_node(parent, child)
                    array(int(node/2)) = parent
                    if (parent == child_1st) then   !parentは入れ替え後なのでもとはchild
                        array(node-1) = child
                    else
                        array(node) = child
                    end if
                end if
            end do

            array_sorted(i) = array(1)
            allocate(temp_array(num_node-i))

            do j = 2, num_node  !入力配列の最上親ノードを抜き取って1つずらした配列作成
                temp_array(j-1) = array(j)
            end do
            
            deallocate(array)
            allocate(array(num_node-i))

            array = temp_array

            deallocate(temp_array)

            num_node_heap = size(array)

        end do

    end subroutine

end module sort_m