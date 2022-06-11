module kdTree_m
    use unstructuredGrid_mod
    use sort_m
    implicit none

    contains

    subroutine cellDivider(grid, after, leftChild, kdTree_cell, rightChild)
        type(unstructuredGrid), intent(in) :: grid
        type(content_t), intent(in) :: after(:)
        type(content_t), intent(out), allocatable :: leftChild(:)
        type(content_t), intent(out), allocatable :: rightChild(:)
        type(content_t), intent(out), allocatable :: kdTree_cell(:)
        integer i, centerID, left_size, right_size

        allocate(kdTree_cell(size(after)))
        
        centerID = int(size(after)/2)+1
        kdTree_cell(1)%originID = after(centerID)%originID
        kdTree_cell(1)%value = grid%CELLs(after(centerID)%originID)%center(1)

        left_size = centerID-1
        allocate(leftChild(left_size))
        do i = 1, left_size
            leftChild(i)%originID = after(i)%originID
            leftChild(i)%value = grid%CELLs(after(i)%originID)%center(1)
        end do

        if(mod(size(after),2) == 0) then
            right_size = centerID-2
            allocate(rightChild(right_size))
            do i = 1, right_size
                rightChild(i)%originID = after(i+centerID)%originID
                rightChild(i)%value = grid%CELLs(after(i+centerID)%originID)%center(1)
            end do   
        else
            right_size = centerID-1
            allocate(rightChild(right_size))
            do i = 1, right_size
                rightChild(i)%originID = after(i+centerID)%originID
                rightChild(i)%value = grid%CELLs(after(i+centerID)%originID)%center(1)
            end do
        end if

    end subroutine

end module