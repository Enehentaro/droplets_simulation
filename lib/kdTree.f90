module kdTree_m
    use sort_m
    implicit none
    private

    type, public :: node_in_kdTree_t
        private
        integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0, cell_ID = 0
        integer depth
        integer, allocatable :: ID_array(:)
    end type

    public create_kdtree, search_kdtree

    contains

    subroutine create_kdtree(xyz_origin, kdTree)
        real, intent(in) :: xyz_origin(:,:) !セル重心座標配列
        type(content_t), allocatable :: x_origin(:), y_origin(:), z_origin(:)
        type(node_in_kdTree_t), intent(out), allocatable :: kdTree(:)
        type(content_t), allocatable :: array_pre(:), array_sorted(:)
        integer, allocatable :: leftChildIDArray(:), rightChildIDArray(:)
        integer centerID, i, num_node
        integer parentID, child1ID, child2ID
        integer depth, ID_counter

        ID_counter = 1

        num_node = size(xyz_origin, dim=2)
        allocate(kdTree(num_node))

        !各軸に対してコンテンツ配列に
        x_origin = real2content(xyz_origin(1,:))
        y_origin = real2content(xyz_origin(2,:))
        z_origin = real2content(xyz_origin(3,:))

        kdTree(1)%ID_array = x_origin(:)%originID
        kdTree(1)%depth = 0 !最初は深さゼロ

        do i = 1, num_node
            parentID = i

            depth = kdTree(i)%depth

            !各軸を切り替えながら、コンテンツ配列から要素を抽出
            select case(mod(depth, 3))
            case(0)
                array_pre = x_origin(kdTree(i)%ID_array)
            case(1)
                array_pre = y_origin(kdTree(i)%ID_array)
            case(2)
                array_pre = z_origin(kdTree(i)%ID_array)
            end select

            array_sorted = array_pre !←配列サイズを揃えるため
            call heap_sort(array_pre, array_sorted)

            centerID = int(size(array_sorted)/2)+1
            kdTree(i)%cell_ID = array_sorted(centerID)%originID     !ヒープソート結果の中央値
            leftChildIDArray = array_sorted(:centerID-1)%originID   !左側配列のIDだけ取り出す
            rightChildIDArray = array_sorted(centerID+1:)%originID  !右側配列のIDだけ取り出す

            if(size(leftChildIDArray) >= 1) then 
                ID_counter = ID_counter + 1
                child1ID = ID_counter
                kdtree(child1ID)%ID_array = leftChildIDArray
                call set_relation(kdTree, parentID, child1ID, 'left')
            end if

            if(size(rightChildIDArray) >= 1) then
                ID_counter = ID_counter + 1
                child2ID = ID_counter
                kdTree(child2ID)%ID_array = rightChildIDArray
                call set_relation(kdTree, parentID, child2ID, 'right')
            end if

        end do

    end subroutine

    subroutine set_relation(array, parent_ID, child_ID, lr)
        type(node_in_kdTree_t) :: array(:) 
        integer, intent(in) :: parent_ID, child_ID
        character(*), intent(in) :: lr

        array(child_ID)%parent_ID = parent_ID
        array(child_ID)%depth = array(parent_ID)%depth + 1

        select case(lr)
        case('left')
            array(parent_ID)%child_ID_1 = child_ID
        case('right')
            array(parent_ID)%child_ID_2 = child_ID
        case default
            print '("relation ERROR")'
            stop
        end select

    end subroutine

    !ただ根ノードから葉ノードまで一方的に下っているだけなので、未完成
    subroutine search_kdtree(xyz, kdTree, droplet_position, nearest_ID)
        real, intent(in) :: xyz(:,:) 
        type(node_in_kdTree_t), intent(in), allocatable :: kdTree(:)
        real, intent(in) :: droplet_position(3)
        integer depth, switch, parentID, nextChildID
        integer, intent(out) :: nearest_ID

        parentID = 1 

        do
            depth = kdTree(parentID)%depth
            switch = mod(depth,3)+1
            if(droplet_position(switch) <= xyz(switch, kdTree(parentID)%cell_ID)) then 
                nextChildID = kdTree(parentID)%child_ID_1
            else 
                nextChildID = kdTree(parentID)%child_ID_2
            end if

            if(nextChildID == 0) then
                exit
            else
                parentID = nextChildID
            end if

        end do

        nearest_ID = kdTree(parentID)%cell_ID

        print*, droplet_position
        print*, xyz(:, nearest_ID)
        print*, 'parentID =', parentID

    end subroutine

end module
