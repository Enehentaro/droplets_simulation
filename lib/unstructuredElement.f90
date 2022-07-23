module unstructuredElement_m
    !!非構造要素モジュール
    !!非構造格子のベースとなる構造体を定義している
    implicit none
    private

    !>節点構造体（ただの座標配列）
    type, public :: node_t
        real coordinate(3)
    end type

    !>セル構造体（節点のID配列）
    type, public :: cell_t
        integer, allocatable :: nodeID(:)
    end type

    public get_MinMaxCDN, get_cellCenters

    contains

    ! function points2nodeArray(points) result(node_array)
    !     real, intent(in) :: points(:,:)
    !     type(node_t), allocatable :: node_array(:)
    !     integer i, num_points

    !     num_points = size(points, dim=2)
    !     allocate(node_array(num_points))

    !     do i = 1, num_points
    !         node_array(i)%coordinate = points(:,i)
    !     end do

    ! end function

    function get_MinMaxCDN(node) result(MinMax)
        !!節点群の座標の最大最小を返す
        class(node_t), intent(in) :: node(:)
        real MinMax(3,2)

        MinMax(1,1) = minval(node(:)%coordinate(1))
        MinMax(2,1) = minval(node(:)%coordinate(2))
        MinMax(3,1) = minval(node(:)%coordinate(3))
        print*, 'MIN_coordinates=', MinMax(:,1)

        MinMax(1,2) = maxval(node(:)%coordinate(1))
        MinMax(2,2) = maxval(node(:)%coordinate(2))
        MinMax(3,2) = maxval(node(:)%coordinate(3))
        print*, 'MAX_coordinates=', MinMax(:,2)

    end function

    function get_cellCenters(node, cell) result(centers)
        !!すべてのセル重心を計算し、配列で返す
        class(node_t), intent(in) :: node(:)
        class(cell_t), intent(in) :: cell(:)
        real, allocatable :: centers(:,:)
        real x(3)
        integer i, k, num_node, num_cell, nodeID

        num_cell = size(cell)
        allocate(centers(3, num_cell))

        do i = 1, num_cell
            x(:) = 0.0
            num_node = size(cell(i)%nodeID)
            do k = 1, num_node
                nodeID = cell(i)%nodeID(k)
                x(:) = x(:) + node(nodeID)%coordinate(:)
            end do
            centers(:,i) = x(:) / real(num_node)
        end do

    end function
    
end module unstructuredElement_m