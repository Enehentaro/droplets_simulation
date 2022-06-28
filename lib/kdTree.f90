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
        real, intent(in) :: xyz_origin(:,:) !セル重心座標配列(3, num_cell)
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
            ! print '(*(i0, x))', array_sorted(:)%originID
            ! print '(*(g0, x))', array_sorted(:)%value

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

        ! call print_tree(kdTree, xyz_origin)
        call saveAsDOT(kdTree, xyz_origin)

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
        type(node_in_kdTree_t), intent(in) :: kdTree(:)
        real, intent(in) :: droplet_position(3)
        integer depth, switch, parentID, nextChildID
        integer, intent(out) :: nearest_ID
        real mindist
        logical, allocatable :: NotYetCompared(:)

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
        
        allocate(NotYetCompared(size(kdTree)))
        NotYetCompared(:) = .true.

        mindist = norm2(xyz(:,kdTree(parentID)%cell_ID)-droplet_position(:))
        NotYetCompared(parentID) = .false.

        ! print*, 'parentID =', parentID
        ! print*, NotYetCompared

        do
            parentID = kdTree(parentID)%parent_ID
            NotYetCompared(parentID) = .false.
            depth = kdTree(parentID)%depth
            switch = mod(depth,3)+1
            ! print*, 'switch = ', switch
            ! print*, 'parentID =', parentID
            ! print*, 'xyz(:,kdTree(parentID)%cell_ID) =', xyz(:,kdTree(parentID)%cell_ID)
            ! print*, 'droplet_position =', droplet_position
            ! print*, 'kdTree(parentID)%child_ID_1 =', kdTree(parentID)%child_ID_1
            ! print*, 'mindist =', mindist
            ! print*, 'abs(xyz(switch,kdTree(parentID)%cell_ID)-droplet_position(switch)) =', &
            !         abs(xyz(switch,kdTree(parentID)%cell_ID)-droplet_position(switch))

            if(mindist <= abs(xyz(switch,kdTree(parentID)%cell_ID)-droplet_position(switch))) then
                if(NotYetCompared(kdTree(parentID)%child_ID_1)) then
                    NotYetCompared(kdTree(parentID)%child_ID_1) = .false.
                else
                    NotYetCompared(kdTree(parentID)%child_ID_2) = .false.
                end if
            else
                if(NotYetCompared(kdTree(parentID)%child_ID_1)) then
                    ! mindist = min(mindist, norm2(xyz(:, kdTree(parentID)%cell_ID)-droplet_position(:)),&
                    ! norm2(xyz(:, kdTree(parentID)%child_ID_1)-droplet_position(:)))
                    ! NotYetCompared(kdTree(parentID)%child_ID_1) = .false.
                else
                    ! mindist = min(mindist, norm2(xyz(:, kdTree(parentID)%cell_ID)-droplet_position(:)),&
                    ! norm2(xyz(:, kdTree(parentID)%child_ID_2)-droplet_position(:)))
                    ! NotYetCompared(kdTree(parentID)%child_ID_2) = .false.
                end if
            end if

            print*, NotYetCompared
            
            stop

        end do

    end subroutine

    ! subroutine print_tree(kdTree, xyz)
    !     type(node_in_kdTree_t), intent(in) :: kdTree(:)
    !     real, intent(in) :: xyz(:,:)
    !     integer i

    !     print '("=======================================================")'
    !     do i = 1, size(kdTree)
    !         print*, 'ID_in_tree:', i
    !         print*, 'parentID:', kdTree(i)%parent_ID
    !         print*, 'childrenID:', kdTree(i)%child_ID_1, kdTree(i)%child_ID_2
    !         print*, 'cell:', kdTree(i)%cell_ID, xyz(:, kdTree(i)%cell_ID)
    !         print '("=======================================================")'
    !     end do

    ! end subroutine

    subroutine saveAsDOT(kdTree, xyz)
        type(node_in_kdTree_t), intent(in) :: kdTree(:)
        real, intent(in) :: xyz(:,:)
        integer n_unit
        integer i
        character(1), parameter :: dq = '"'

        open(newunit=n_unit, file='Test_check/kdTree.dot')
        write(n_unit, '(A)') 'graph {'

        write(n_unit, '(4x, A)') 'node ['
        write(n_unit, '(2(4x), A)') 'shape = record,'
        write(n_unit, '(4x, A)') '];'
          
        write(n_unit, '()')

        ! node define
        do i = 1, size(kdTree)
            write(n_unit, '(4x, i0, "[label = ", A, "{", i0, "| cell ID : ", i0, "|", 3(f10.5), "}", A, "];")') &
                i, dq, i, kdTree(i)%cell_ID, xyz(:, kdTree(i)%cell_ID), dq
        end do

        write(n_unit, '()')

        ! edge define
        do i = 1, size(kdTree)
            if(kdTree(i)%child_ID_1 /= 0) then
                write(n_unit, '(4x, i0, " -- ", i0, ";")') i, kdTree(i)%child_ID_1
            end if
            if(kdTree(i)%child_ID_2 /= 0) then
                write(n_unit, '(4x, i0, " -- ", i0, ";")') i, kdTree(i)%child_ID_2
            end if
        end do


        write(n_unit, '(A)') '}'

    end subroutine

end module
