module kdTree_m
    use unstructuredGrid_mod
    use sort_m
    implicit none

    type nodeTree_t
        integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0 ,cell_ID = 0
        type(content_t), allocatable :: left_array(:), right_array(:)
    end type

    ! integer :: ID_counter = 1

    contains

    subroutine create_kdtree(before)
        type(content_t), intent(in) :: before(:)
        type(content_t), allocatable :: after(:)
        type(nodeTree_t), allocatable :: node_tree(:)
        type(content_t), allocatable :: pre_leftChild(:), pre_rightChild(:), leftChild(:), rightChild(:)
        integer centerID, switch, left_size, right_size, i, arraySize, n, num_loop
        integer child1ID, child2ID
        integer :: depth = 1
        integer :: parentID = 1
        logical :: left_flag = .true.

        n = int(log10(dble(size(before)+1)) / log10(2.0))
        num_loop = 2**n - 1
        allocate(node_tree(2*num_loop+1))

        switch = mod(depth-1,3)+1
        call heap_sort(before, switch, after)

        centerID = int(size(after)/2)+1
        node_tree(1)%cell_ID = after(centerID)%originID ! ヒープソートの1回目の結果の中央値 

        call split_cell(after(:), pre_leftChild, pre_rightChild)
        depth = depth + 1

        switch = mod(depth-1,3)+1

        if(size(pre_leftChild) >= 1) then 
            call heap_sort(pre_leftChild, switch, leftChild)
            deallocate(pre_leftChild)
            child1ID = 2*parentID
            allocate(node_tree(1)%left_array(size(leftChild)))
            node_tree(1)%left_array(:) = leftChild(:)
        end if

        if(size(pre_rightChild) >= 1) then
            call heap_sort(pre_rightChild, switch, rightchild)
            deallocate(pre_rightChild)
            child2ID = 2*parentID + 1
            allocate(node_tree(1)%right_array(size(rightChild)))
            node_tree(1)%right_array(:) = rightChild(:)
        end if

        centerID = int(size(leftChild)/2)+1
        node_tree(child1ID)%cell_ID = leftChild(centerID)%originID
        deallocate(leftChild)
        centerID = int(size(rightChild)/2)+1
        node_tree(child2ID)%cell_ID = rightChild(centerID)%originID
        deallocate(rightChild)

        call solve_relation(node_tree,parentID,child1ID,child2ID)

        depth = 3

        ! deallocate(before)
        deallocate(after)
        print*, parentID
        print*, size(node_tree(parentID)%left_array)
        print*, size(node_tree(parentID)%right_array)

        do i = 2, num_loop
            parentID = i
            
            if(left_flag) then

                print*, "parentID = ", parentID

                left_flag = .false.
                if(size(node_tree(i/2)%left_array) /= 1) then

                    call split_cell(node_tree(i/2)%left_array, pre_leftChild, pre_rightChild)

                    switch = mod(depth-1,3)+1

                    call heap_sort(pre_leftChild, switch, leftChild)
                    deallocate(pre_leftChild)
                    child1ID = 2*parentID
                    allocate(node_tree(parentID)%left_array(size(leftChild)))
                    node_tree(parentID)%left_array(:) = leftChild(:)
                    print*, "child1ID = ", child1ID

                    call heap_sort(pre_rightChild, switch, rightchild)
                    deallocate(pre_rightChild)
                    child2ID = 2*parentID + 1
                    allocate(node_tree(parentID)%right_array(size(rightChild)))
                    node_tree(parentID)%right_array(:) = rightChild(:)
                    print*, "child2ID = ", child2ID
                    
                    centerID = int(size(leftChild)/2)+1
                    node_tree(child1ID)%cell_ID = leftChild(centerID)%originID
                    deallocate(leftChild)
                    centerID = int(size(rightChild)/2)+1
                    node_tree(child2ID)%cell_ID = rightChild(centerID)%originID
                    deallocate(rightChild)

                    call solve_relation(node_tree,parentID,child1ID,child2ID)

                    print*, "parentID = ", parentID
                    print*, size(node_tree(parentID)%left_array)
                    print*, size(node_tree(parentID)%right_array)

                end if
                
            else

                print*, "parentID = ", parentID
                left_flag = .true.
                if(size(node_tree(int(i/2))%right_array) /= 1) then

                    call split_cell(node_tree(int(i/2))%right_array, pre_leftChild, pre_rightChild)
                    switch = mod(depth-1,3)+1

                    call heap_sort(pre_leftChild, switch, leftChild)
                    deallocate(pre_leftChild)
                    child1ID = 2*parentID
                    allocate(node_tree(parentID)%left_array(size(leftChild)))
                    node_tree(parentID)%left_array(:) = leftChild(:)
                    print*, "child1ID = ", child1ID

                    call heap_sort(pre_rightChild, switch, rightchild)
                    deallocate(pre_rightChild)
                    child2ID = 2*parentID + 1
                    allocate(node_tree(parentID)%right_array(size(rightChild)))
                    node_tree(parentID)%right_array(:) = rightChild(:)
                    print*, "child2ID = ", child2ID
                    
                    centerID = int(size(leftChild)/2)+1
                    node_tree(child1ID)%cell_ID = leftChild(centerID)%originID
                    deallocate(leftChild)
                    centerID = int(size(rightChild)/2)+1
                    node_tree(child2ID)%cell_ID = rightChild(centerID)%originID
                    deallocate(rightChild)

                    call solve_relation(node_tree,parentID,child1ID,child2ID)

                    print*, "parentID = ", parentID
                    print*, size(node_tree(parentID)%left_array)
                    print*, size(node_tree(parentID)%right_array)

                    deallocate(node_tree(int(i/2))%left_array)
                    deallocate(node_tree(int(i/2))%right_array)

                end if

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
        type(nodeTree_t):: array(:) 
        integer :: parent_ID,child_ID_1,child_ID_2

        array(parent_ID)%child_ID_1 = child_ID_1 
        array(parent_ID)%child_ID_2 = child_ID_2 

        array(child_ID_1)%parent_ID = parent_ID 
        array(child_ID_2)%parent_ID = parent_ID 

    end subroutine

end module
