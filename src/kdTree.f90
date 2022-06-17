module kdTree_m
    ! use unstructuredGrid_mod
    ! use sort_m
    implicit none
    private

    ! type node
    !     integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0 ,cell_ID = 0 
    ! end type 

    ! integer ID_counter = 1

    ! contains

    ! call heap_sort()
    ! node_tree(1)%cell_ID = !ヒープソートの1回目の結果の中央値    

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


    ! subroutine solve_relation(array,parent_ID,child_ID_1,child_ID_2)
    !     type(nonde):: array(:) 
    !     integer :: parent_ID,child_ID_1,child_ID_2

    !     array(parent_ID)%child_ID_1 = child_ID_1 
    !     array(parent_ID)%child_ID_2 = child_ID_2 

    !     array(child_ID_1)%parent_ID = parent_ID 
    !     array(child_ID_2)%parent_ID = parent_ID 

    ! end subroutine solve_relation

end module