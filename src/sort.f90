module sort_m
    implicit none

    type content
        integer cellID
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
    
    subroutine heap_sort(array_origin, array_sorted)
        type(content), intent(in) :: array_origin(:)
        type(content), intent(out), allocatable :: array_sorted(:)

        type(content), allocatable :: calc_array(:)
        type(content), allocatable :: temp_array(:)  !配列の要素を減らすための一時的な配列
        type(content) parent, child, larger_child, child_1st, child_2nd
        integer i, j, node, num_node_heap
        integer parentID, largerChildID

        num_node_heap = size(array_origin)   !初めは要素全体のノード数とヒープ構造のノード数が同じ
        allocate(calc_array(size(array_origin)))
        allocate(array_sorted(size(array_origin)))
        
        do i = 1, size(array_origin)
            calc_array(i)%cellID = array_origin(i)%cellID
            calc_array(i)%axis = array_origin(i)%axis
        end do

        do i = 1, size(array_origin)  !ソート後の配列に格納するループ

            if (mod(num_node_heap,2) == 0) then  !要素数の偶奇判定(偶数のときだけ特別な処理)
                parent = calc_array(num_node_heap/2)
                child = calc_array(num_node_heap)

                if(parent%axis < child%axis) then !要素数が偶数のとき末端のノードだけ2分木にならないのでその処理
                    call swap_content(calc_array, num_node_heap/2, num_node_heap)
                end if

                num_node_heap = num_node_heap-1   !末端の1分木の分だけ比較に使用するヒープ構造の要素数を減らす
            
            end if
                
            do node = num_node_heap, 3, -2   !ヒープソート一回分のループ
                parentID = int(node/2)
                largerChildID = get_largerID(calc_array(:)%axis, node-1, node)

                if(calc_array(parentID)%axis < calc_array(largerChildID)%axis) then
                    call swap_content(calc_array, parentID, largerChildID)
                end if
            end do

            array_sorted(i) = calc_array(1)
            allocate(temp_array(size(calc_array)-1))

            do j = 2, size(calc_array)  !入力配列の最上親ノードを抜き取って1つずらした配列作成
                temp_array(j-1) = calc_array(j)
            end do

            deallocate(calc_array)
            allocate(calc_array(size(array_origin)-i))

            calc_array = temp_array

            deallocate(temp_array)

            num_node_heap = size(calc_array)

        end do

    end subroutine

end module sort_m