module kdTree_m
    use unstructuredGrid_mod
    use sort_m
    implicit none

    type nodeTree_t
        integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0 ,cell_ID = 0 
    end type

    integer :: ID_counter = 1

    contains

    subroutine create_kdtree(before)
        type(content_t), intent(in) :: before(:)
        type(nodeTree_t), allocatable :: node_tree(:)
        type(content_t), allocatable :: pre_leftChild(:), pre_rightChild(:), leftChild(:), rightChild(:)
        integer centerID, left_size, right_size

        call kdtree_setup(before, leftChild, rightChild)    

        ! do i = 1, num_ID-1
        !     parent_ID = i

        !     pre_leftChild(:) = 
        !     pre_rightChild(:) = 

        !     if(size(pre_leftChild) >= 1) then 
        !         call heap_sort(pre_leftChild, leftchild)
        !         ID_counter = ID_counter + 1
        !         child_ID_1 = ID_counter
        !     end if

        !     if(size(pre_rightChild) >= 1) then
        !         call heap_sort(pre_rightChild, rightchild)
        !         ID_counter = ID_counter + 1
        !         child_ID_2 = ID_counter
        !     end if

        !     node_tree(child_ID_1)%cell_ID = !上記の中央値
        !     node_tree(child_ID_2)%cell_ID = !上記の中央値
        !     call solve_relation(node_tree,parent_ID,child_ID_1,child_ID_2)

        ! end do

    end subroutine

    subroutine kdtree_setup(before, leftChild, rightChild)
        type(content_t), intent(in) :: before(:)
        type(content_t), intent(out), allocatable :: leftChild(:), rightChild(:)
        type(content_t), allocatable :: after(:)
        type(nodeTree_t), allocatable :: node_tree(:)
        type(content_t), allocatable :: pre_leftChild(:), pre_rightChild(:)
        integer centerID, child_centerID, left_size, right_size, child_ID_1, child_ID_2
        integer i, n_unit
        integer :: parent_ID = 1
        integer :: depth = 1

        allocate(after(size(before)))
        call heap_sort(before(:), depth, after(:))
        
        allocate(node_tree(size(after)))
        centerID = int(size(after)/2)+1

        node_tree(1)%cell_ID = after(centerID)%originID ! ヒープソートの1回目の結果の中央値
        
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

        open(newunit = n_unit, file = "Test/first_divide.txt", status = 'replace')
            do i = 1, size(pre_leftChild) 
                write(n_unit,'(I0)', advance='no') pre_leftChild(i)%originID
                write(n_unit,'(3(f12.5))') pre_leftChild(i)%coordinate(1), pre_leftChild(i)%coordinate(2), pre_leftChild(i)%coordinate(3)
            end do
            write(n_unit,'(A)')
            write(n_unit,'(I0)') node_tree(1)%cell_ID
            write(n_unit,'(A)')
            do i = 1, size(pre_rightChild)
                write(n_unit,'(I0)', advance='no') pre_rightChild(i)%originID
                write(n_unit,'(3(f12.5))') pre_rightChild(i)%coordinate(1), pre_rightChild(i)%coordinate(2), pre_rightChild(i)%coordinate(3)
            end do
        close(n_unit)

        depth = depth + 1

        if(size(pre_leftChild) >= 1) then
            allocate(leftChild(size(pre_leftChild)))
            call heap_sort(pre_leftChild, depth, leftchild)
            ID_counter = ID_counter + 1
            child_ID_1 = ID_counter
        end if
        if(size(pre_rightChild) >= 1) then
            allocate(rightChild(size(pre_rightChild)))
            call heap_sort(pre_rightChild, depth, rightchild)
            ID_counter = ID_counter + 1
            child_ID_2 = ID_counter
        end if

        child_centerID = int(size(leftChild)/2)+1
        node_tree(child_ID_1)%cell_ID = leftchild(child_centerID)%originID ! 左子ノードの中央値
        node_tree(child_ID_2)%cell_ID = rightchild(child_centerID)%originID ! 右子ノードの中央値
        call solve_relation(node_tree,parent_ID,child_ID_1,child_ID_2)

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
