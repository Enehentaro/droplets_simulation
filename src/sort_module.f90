module sort_mod
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
    

    subroutine construct_heap(array, array_sorted)
        real, intent(in), allocatable :: array(:)
        real, intent(out), allocatable :: array_sorted(:)
        real, allocatable :: temp_array(:)  !配列の要素を減らすための一時的な配列
        real parent, child, child_1st, child_2nd
        integer i, j, node, num_node

        num_node = size(array)
        allocate(array_sorted(numnode))

        if (mod(num_node,2) == 0) then  !要素数の偶奇判定
            parent = array(num_node/2)
            child = array(num_node)
            if(parent < child) then
                call swap_node(parent, child)
                array(num_node/2) = parent
                array(num_node) = child
            end if
        end if

        do i = 1, num_node  !ソート後の配列に格納するループ
            do node = num_node, 3, -2   !ヒープソート一回分のループ
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

        end do
    end subroutine
        

    subroutine heap_sort
    end subroutine

end module sort_mod