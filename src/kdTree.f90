module kdTree_m
    use unstructuredGrid_mod
    use sort_m
    implicit none

    contains

    subroutine cellDivider(grid, after, leftChild, kdTree_cell, rightChild)
        type(unstructuredGrid), intent(in) :: grid
        type(content), intent(in) :: after(:)
        type(content), intent(out), allocatable :: leftChild(:)
        type(content), intent(out), allocatable :: rightChild(:)
        type(content), intent(out), allocatable :: kdTree_cell(:)
        integer i, centerID, left_size, right_size

        allocate(kdTree_cell(size(after)))
        
        centerID = int(size(after)/2)+1
        kdTree_cell(1)%cellID = after(centerID)%cellID
        kdTree_cell(1)%axis = grid%CELLs(after(centerID)%cellID)%center(1)

        left_size = centerID-1
        allocate(leftChild(left_size))
        do i = 1, left_size
            leftChild(i)%cellID = after(i)%cellID
            leftChild(i)%axis = grid%CELLs(after(i)%cellID)%center(1)
        end do

        if(mod(size(after),2) == 0) then
            right_size = centerID-2
            allocate(rightChild(right_size))
            do i = 1, right_size
                rightChild(i)%cellID = after(i+centerID)%cellID
                rightChild(i)%axis = grid%CELLs(after(i+centerID)%cellID)%center(1)
            end do   
        else
            right_size = centerID-1
            allocate(rightChild(right_size))
            do i = 1, right_size
                rightChild(i)%cellID = after(i+centerID)%cellID
                rightChild(i)%axis = grid%CELLs(after(i+centerID)%cellID)%center(1)
            end do
        end if

    end subroutine

end module