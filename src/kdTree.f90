module kdTree_m
    ! use unstructuredGrid_mod
    ! use sort_m
    implicit none
    private

    type nodeTree_t
        integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0 ,cell_ID = 0
        type(content_t), allocatable :: array(:)
    end type

    integer :: ID_counter = 1

    contains

    subroutine create_kdtree(before, node_tree)
        type(content_t), intent(in) :: before(:)
        type(content_t), allocatable :: after(:)
        type(nodeTree_t), intent(out), allocatable :: node_tree(:)
        type(content_t), allocatable :: pre_leftChild(:), pre_rightChild(:), leftChild(:), rightChild(:)
        integer centerID, switch, left_size, right_size, i, arraySize
        integer parentID, child1ID, child2ID
        integer depth

        depth = 1
        parentID = 1

        arraySize = size(before)
        allocate(node_tree(arraySize))

        switch = mod(depth-1,3)+1
        call heap_sort(before, switch, after)

        node_tree(1)%array = after(:)

        centerID = int(size(after)/2)+1
        node_tree(1)%cell_ID = after(centerID)%originID ! ヒープソートの1回目の結果の中央値 

        call split_cell(after(:), pre_leftChild, pre_rightChild)
        depth = depth + 1

        switch = mod(depth-1,3)+1

        if(size(pre_leftChild) >= 1) then 
            call heap_sort(pre_leftChild, switch, leftChild)
            deallocate(pre_leftChild)
            ID_counter = ID_counter + 1
            child1ID = ID_counter
            allocate(node_tree(child1ID)%array(size(leftChild)))
            node_tree(child1ID)%array(:) = leftChild(:)
        end if

        if(size(pre_rightChild) >= 1) then
            call heap_sort(pre_rightChild, switch, rightchild)
            deallocate(pre_rightChild)
            ID_counter = ID_counter + 1
            child2ID = ID_counter
            allocate(node_tree(child2ID)%array(size(rightChild)))
            node_tree(child2ID)%array(:) = rightChild(:)
        end if

        centerID = int(size(leftChild)/2)+1
        node_tree(child1ID)%cell_ID = leftChild(centerID)%originID
        deallocate(leftChild)
        centerID = int(size(rightChild)/2)+1
        node_tree(child2ID)%cell_ID = rightChild(centerID)%originID
        deallocate(rightChild)

        call solve_relation(node_tree,parentID,child1ID,child2ID)

        depth = 3

        deallocate(after)

        do i = 2, arraySize
            parentID = i

            if(size(node_tree(parentID)%array) /= 1) then

                call split_cell(node_tree(parentID)%array, pre_leftChild, pre_rightChild)

                switch = mod(depth-1,3)+1

                if(size(pre_leftChild) >= 1) then 
                    call heap_sort(pre_leftChild, switch, leftChild)
                    deallocate(pre_leftChild)
                    ID_counter = ID_counter + 1
                    child1ID = ID_counter
                    allocate(node_tree(child1ID)%array(size(leftChild)))
                    node_tree(child1ID)%array(:) = leftChild(:)

                    centerID = int(size(leftChild)/2)+1
                    node_tree(child1ID)%cell_ID = leftChild(centerID)%originID
                    deallocate(leftChild)
                end if

                if(size(pre_rightChild) >= 1) then
                    call heap_sort(pre_rightChild, switch, rightchild)
                    deallocate(pre_rightChild)
                    ID_counter = ID_counter + 1
                    child2ID = ID_counter
                    allocate(node_tree(child2ID)%array(size(rightChild)))
                    node_tree(child2ID)%array(:) = rightChild(:)
                    
                    centerID = int(size(rightChild)/2)+1
                    node_tree(child2ID)%cell_ID = rightChild(centerID)%originID
                    deallocate(rightChild)
                else
                    child2ID = child1ID
                end if
            
                call solve_relation(node_tree,parentID,child1ID,child2ID)

                deallocate(node_tree(parentID)%array)

            end if

            if(i == 2**(depth-1)-1) then
                depth = depth + 1
            end if

        end do

    end subroutine

    subroutine split_cell(after, pre_leftChild, pre_rightChild)
        type(content_t), intent(in) :: after(:)
        type(content_t), intent(out), allocatable :: pre_leftChild(:), pre_rightChild(:)
        integer left_size, right_size, i, centerID

        centerID = int(size(after)/2)+1
        left_size = centerID-1
        allocate(pre_leftChild(left_size))
        do i = 1, left_size
            pre_leftChild(i) = after(i)
        end do

        ! 右子ノードサイズの偶奇で場合分け
        if(mod(size(after),2) == 0) then
            right_size = centerID-2
            allocate(pre_rightChild(right_size))
            do i = 1, right_size
                pre_rightChild(i) = after(i+centerID)
            end do
        else
            right_size = centerID-1
            allocate(pre_rightChild(right_size))
            do i = 1, right_size
                pre_rightChild(i) = after(i+centerID)
            end do
        end if

    end subroutine

    subroutine solve_relation(array,parent_ID,child_ID_1,child_ID_2)
        type(nodeTree_t) :: array(:) 
        integer, intent(in) :: parent_ID,child_ID_1,child_ID_2

        if(child_ID_1 == child_ID_2) then
            array(parent_ID)%child_ID_1 = child_ID_1
            array(parent_ID)%child_ID_2 = child_ID_1
            array(child_ID_1)%parent_ID = parent_ID
        else
            array(parent_ID)%child_ID_1 = child_ID_1 
            array(parent_ID)%child_ID_2 = child_ID_2 

            array(child_ID_1)%parent_ID = parent_ID 
            array(child_ID_2)%parent_ID = parent_ID 
        end if
    end subroutine

    subroutine search_kdtree(before,node_tree, droplet_position, nearest_ID)
        type(content_t),intent(in),allocatable :: before(:) 
        type(nodeTree_t), intent(in), allocatable :: node_tree(:)
        real, intent(in) :: droplet_position(3)
        integer nearest_ID ,depth ,switch ,n ,parentID 

        parentID = 1 
        depth = 1 

        do while (size(node_tree(parentID)%array) > 1)
            switch = mod(depth-1,3)+1
            depth = depth + 1
            if(droplet_position(switch) <= before(node_tree(parentID)%cell_ID)%coordinate(switch)) then 
                parentID = node_tree(parentID)%child_ID_1
            else 
                parentID = node_tree(parentID)%child_ID_2
            end if
        end do
        print*, 'parentID =', parentID 

    end subroutine

end module
