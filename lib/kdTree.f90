module kdTree_m
    use sort_m
    implicit none
    private

    type node_in_kdTree_t
        private
        integer :: parent_ID = 0, child_ID_1 = 0, child_ID_2 = 0, cell_ID = 0
        integer depth
        integer, allocatable :: ID_array(:)
    end type

    type, public :: kdTree
        private
        type(node_in_kdTree_t), allocatable :: node(:)
        contains
        procedure set_relation, saveAsDOT, create_childlist
        procedure :: search => search_kdtree
    end type

    public kdTree_

    contains

    type(kdTree) function kdTree_(xyz_origin)
        real, intent(in) :: xyz_origin(:,:) !セル重心座標配列(3, num_cell)
        type(content_t), allocatable :: x_origin(:), y_origin(:), z_origin(:)
        type(content_t), allocatable :: array_pre(:), array_sorted(:)
        integer, allocatable :: leftChildIDArray(:), rightChildIDArray(:)
        integer centerID, i, num_node
        integer parentID, child1ID, child2ID
        integer depth, ID_counter

        ID_counter = 1

        num_node = size(xyz_origin, dim=2)
        allocate(kdTree_%node(num_node))

        !各軸に対してコンテンツ配列に
        x_origin = real2content(xyz_origin(1,:))
        y_origin = real2content(xyz_origin(2,:))
        z_origin = real2content(xyz_origin(3,:))

        kdTree_%node(1)%ID_array = x_origin(:)%originID
        kdTree_%node(1)%depth = 0 !最初は深さゼロ

        do i = 1, num_node
            parentID = i

            depth = kdTree_%node(i)%depth

            !各軸を切り替えながら、コンテンツ配列から要素を抽出
            select case(mod(depth, 3))
            case(0)
                array_pre = x_origin(kdTree_%node(i)%ID_array)
            case(1)
                array_pre = y_origin(kdTree_%node(i)%ID_array)
            case(2)
                array_pre = z_origin(kdTree_%node(i)%ID_array)
            end select

            array_sorted = array_pre !←配列サイズを揃えるため
            call heap_sort(array_pre, array_sorted)
            ! print '(*(i0, x))', array_sorted(:)%originID
            ! print '(*(g0, x))', array_sorted(:)%value

            centerID = int(size(array_sorted)/2)+1
            kdTree_%node(i)%cell_ID = array_sorted(centerID)%originID     !ヒープソート結果の中央値
            leftChildIDArray = array_sorted(:centerID-1)%originID   !左側配列のIDだけ取り出す
            rightChildIDArray = array_sorted(centerID+1:)%originID  !右側配列のIDだけ取り出す

            if(size(leftChildIDArray) >= 1) then 
                ID_counter = ID_counter + 1
                child1ID = ID_counter
                kdTree_%node(child1ID)%ID_array = leftChildIDArray
                call kdTree_%set_relation(parentID, child1ID, 'left')
            end if

            if(size(rightChildIDArray) >= 1) then
                ID_counter = ID_counter + 1
                child2ID = ID_counter
                kdTree_%node(child2ID)%ID_array = rightChildIDArray
                call kdTree_%set_relation(parentID, child2ID, 'right')
            end if

        end do

        ! call print_tree(kdTree, xyz_origin)

    end function

    subroutine set_relation(self, parent_ID, child_ID, lr)
        class(kdTree) self
        integer, intent(in) :: parent_ID, child_ID
        character(*), intent(in) :: lr

        self%node(child_ID)%parent_ID = parent_ID
        self%node(child_ID)%depth = self%node(parent_ID)%depth + 1

        select case(lr)
        case('left')
            self%node(parent_ID)%child_ID_1 = child_ID
        case('right')
            self%node(parent_ID)%child_ID_2 = child_ID
        case default
            print '("relation ERROR")'
            stop
        end select

    end subroutine

    !ただ根ノードから葉ノードまで一方的に下っているだけなので、未完成
    subroutine search_kdTree(self, xyz, droplet_position, nearest_ID)
        class(kdTree), intent(in) :: self
        real, intent(in) :: xyz(:,:) 
        real, intent(in) :: droplet_position(3)
        integer depth, switch, parentID, nextChildID, i
        integer, intent(out) :: nearest_ID
        real mindist
        logical, allocatable :: NotYetCompared(:)
        integer, allocatable :: childIDlist(:)

        parentID = 1

        do
            depth = self%node(parentID)%depth
            switch = mod(depth,3)+1
            if(droplet_position(switch) <= xyz(switch, self%node(parentID)%cell_ID)) then 
                nextChildID = self%node(parentID)%child_ID_1
            else 
                nextChildID = self%node(parentID)%child_ID_2
            end if

            if(nextChildID == 0) then
                exit
            else
                parentID = nextChildID
            end if

        end do

        nearest_ID = self%node(parentID)%cell_ID
        
        allocate(NotYetCompared(size(self%node)))
        NotYetCompared(:) = .true.

        mindist = norm2(xyz(:,nearest_ID)-droplet_position(:))
        NotYetCompared(parentID) = .false.

        print*, 'mindist =', mindist

        do
            parentID = self%node(parentID)%parent_ID
            NotYetCompared(parentID) = .false.
            depth = self%node(parentID)%depth
            switch = mod(depth,3)+1

            if(mindist <= abs(xyz(switch,self%node(parentID)%cell_ID)-droplet_position(switch))) then
                if(NotYetCompared(self%node(parentID)%child_ID_1)) then
                    print*, 'before_parentID =', parentID
                    call self%create_childlist(parentID, 'left', childIDlist)
                    do i = 1, size(childIDlist)
                        NotYetCompared(childIDlist(i)) = .false.
                    end do
                else
                    call self%create_childlist(parentID, 'right', childIDlist)
                    print*, 'before_parentID =', parentID
                    do i = 1, size(childIDlist)
                        NotYetCompared(childIDlist(i)) = .false.
                    end do
                end if
            else
                if(NotYetCompared(self%node(parentID)%child_ID_1)) then

                else

                end if
            end if

            print*, NotYetCompared

            print*, "after_parentID =", parentID

            if(parentID == 1) stop

        end do

    end subroutine

    subroutine create_childlist(self, parentID, lr, childIDlist)
        class(kdTree) self
        integer, intent(in) :: parentID
        character(*), intent(in) :: lr
        integer, allocatable, intent(out) :: childIDlist(:)
        integer, allocatable :: tmp_childIDlist(:)
        integer cnt, storage_cnt, origin_parentID, leftID, rightID, i, tmp_parentID

        tmp_parentID = parentID
        ! tmp_childIDlist:親IDと親IDに付随する子供ID全てを一時的に格納する配列
        ! とりあえずtmp_childIDlistの配列の要素数をノード数と合わせる 
        allocate(tmp_childIDlist(size(self%node)))
        tmp_childIDlist = 0
        tmp_childIDlist(1) = tmp_parentID

        ! NotYetComparedがtrueの方の子供IDをtmp_childIDlist(2)に格納
        select case(lr)
            case('left')
                if(self%node(parentID)%child_ID_1 /= 0) then
                    tmp_childIDlist(2) = self%node(tmp_parentID)%child_ID_1
                    tmp_parentID = tmp_childIDlist(2)
                end if
            case('right')
                if(self%node(tmp_parentID)%child_ID_2 /= 0) then
                    tmp_childIDlist(2) = self%node(tmp_parentID)%child_ID_2
                    tmp_parentID = tmp_childIDlist(2)
                end if
        end select

        ! tmp_childIDlist(2)に付随する子供ID全てを格納(左の子供ID＝0 or 左の子供ID＝0の場合も含む)
        ! 左の子供ID＝右の子供ID＝0となったタイミングでループ脱出
        cnt = 2
        do
            tmp_childIDlist(2*cnt-1) = self%node(tmp_childIDlist(cnt))%child_ID_1
            tmp_childIDlist(2*cnt) = self%node(tmp_childIDlist(cnt))%child_ID_2
            if(tmp_childIDlist(2*cnt-1) == 0 .and. tmp_childIDlist(2*cnt) == 0) exit
            cnt = cnt + 1           
        end do

        ! 子供ID=0を含まないchildIDlistのサイズを決定
        cnt = 0
        do i = 1, size(tmp_childIDlist)
            if(tmp_childIDlist(i) == 0 .and. tmp_childIDlist(i+1) == 0) exit
            if(tmp_childIDlist(i) /= 0) then
                cnt = cnt + 1
            end if
        end do

        allocate(childIDlist(cnt))

        ! childIDlistに左詰めで0以外の子供IDを格納していく
        storage_cnt = 1
        do i = 1, cnt
            if(tmp_childIDlist(i) == 0 .and. tmp_childIDlist(i+1) == 0) exit
            if(tmp_childIDlist(i) /= 0) then
                childIDlist(storage_cnt) = tmp_childIDlist(i)
                storage_cnt = storage_cnt + 1
            end if
        end do

        deallocate(tmp_childIDlist)
        
        print*, "childIDlist", childIDlist

    end subroutine

    subroutine saveAsDOT(self, xyz, fname)
        class(kdTree), intent(in) :: self
        real, intent(in) :: xyz(:,:)
        character(*), intent(in) :: fname
        integer n_unit
        integer i
        character(1), parameter :: dq = '"'

        open(newunit=n_unit, file=fname)
        write(n_unit, '("graph {")')

        write(n_unit, '(4x, "node [")')
        write(n_unit, '(2(4x), "shape = record,")')
        write(n_unit, '(4x, "];")')
          
        write(n_unit, '()')

        ! node define
        do i = 1, size(self%node)
            write(n_unit, '(4x, i0, "[label = ", A, "{", i0, "| cell ID : ", i0, "|", 3(f10.5), "}", A, "];")') &
                i, dq, i, self%node(i)%cell_ID, xyz(:, self%node(i)%cell_ID), dq
        end do

        write(n_unit, '()')

        ! edge define
        do i = 1, size(self%node)
            if(self%node(i)%child_ID_1 /= 0) then
                write(n_unit, '(4x, i0, " -- ", i0, ";")') i, self%node(i)%child_ID_1
            end if
            if(self%node(i)%child_ID_2 /= 0) then
                write(n_unit, '(4x, i0, " -- ", i0, ";")') i, self%node(i)%child_ID_2
            end if
        end do


        write(n_unit, '("}")')

    end subroutine

end module
