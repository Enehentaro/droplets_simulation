module CUBE_mod
    use plot3d_operator
    implicit none

    type(cube_inP3D), allocatable, private :: mesh_full(:)
    type(cube_inP3D), allocatable, private :: mesh(:)

    contains

    subroutine read_CUBE_data(fName, path)
        character(*), intent(in) :: fName, path

        if(.not.allocated(mesh_full)) mesh_full = read_plot3d_multigrid(path//'mesh.g')   !Gファイル読み込み
        call read_plot3d_f(mesh_full, fName)   !Fファイル読み込み
        mesh = change_resolution(mesh_full, 4)   !低解像度化

    end subroutine read_CUBE_data

    function get_minMax_CUBE() result(minMax)
        real minMax(6)

        minMax(:) =  minMax_coordinates(mesh)

    end function get_minMax_CUBE

    integer function get_cube_contains(cdn)
        real, intent(in) :: cdn(3)
        integer i_cube, num_cube

        get_cube_contains = 0
        num_cube = size(mesh)

        do i_cube = 1, num_cube
            if  ((cdn(1) >= minval(mesh(i_cube)%nodes(:,:,:)%x)).and. &
                (cdn(2) >= minval(mesh(i_cube)%nodes(:,:,:)%y)).and. &
                (cdn(3) >= minval(mesh(i_cube)%nodes(:,:,:)%z))) then
                    if  ((cdn(1) <= maxval(mesh(i_cube)%nodes(:,:,:)%x)).and. &
                        (cdn(2) <= maxval(mesh(i_cube)%nodes(:,:,:)%y)).and. &
                        (cdn(3) <= maxval(mesh(i_cube)%nodes(:,:,:)%z))) then

                        get_cube_contains = i_cube
                        return

                    end if

            end if

        end do

        print*, 'The Point is Out of Area.'

    end function get_cube_contains

    function nearest_node(cdn, cubeID) result(rnodeID)
        real, intent(in) :: cdn(3)
        integer, intent(in) :: cubeID
        integer rnodeID(3)
        integer i,j,k, i_max,j_max,k_max
        real position(3)
        real, allocatable :: distance(:,:,:)
        type(cube_inP3D) cube_

        cube_ = mesh(cubeID)

        i_max = size(cube_%nodes(:,:,:), dim=1)
        j_max = size(cube_%nodes(:,:,:), dim=2)
        k_max = size(cube_%nodes(:,:,:), dim=3)
        allocate(distance(i_max, j_max, k_max), source=1.e9)
        !$OMP parallel private(position)
        !$OMP do collapse(3)
        do k = 1, k_max
            do j = 1, j_max
                do i = 1, i_max
                    position(1) = cdn(1) - cube_%nodes(i,j,k)%x
                    position(2) = cdn(2) - cube_%nodes(i,j,k)%y
                    position(3) = cdn(3) - cube_%nodes(i,j,k)%z
                    distance(i,j,k) = norm2(position(:))
                end do
            end do
        end do
        !$OMP end do
        !$OMP end parallel

        rnodeID = minloc(distance)

    end function nearest_node

    function nearer_node(cdn, inodeID, cubeID) result(rnodeID)
        real, intent(in) :: cdn(3)
        integer, intent(in) :: inodeID(3)
        integer, intent(in) :: cubeID
        integer nodeID(3), rnodeID(3), i, maxID
        real distance_min, position(3)
        integer,parameter :: diff(3,6) &
                            = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], shape(diff))
        logical update
        type(cube_inP3D) cube_
        type(node_inCUBE) node_

        cube_ = mesh(cubeID)

        node_ = cube_%nodes(inodeID(1), inodeID(2), inodeID(3))
        position(1) = cdn(1) - node_%x
        position(2) = cdn(2) - node_%y
        position(3) = cdn(3) - node_%z
        distance_min = norm2(position(:))
        rnodeID = inodeID
        update = .true.

        maxID = size(cube_%nodes(:,:,:), dim=1)

        do while(update)
            update = .false.

            check : do i = 1, 6

                nodeID(:) = rnodeID(:) + diff(:,i)

                if((minval(nodeID) < 1).or.(maxval(nodeID) > maxID)) cycle check

                node_ = cube_%nodes(nodeID(1), nodeID(2), nodeID(3))
                position(1) = cdn(1) - node_%x
                position(2) = cdn(2) - node_%y
                position(3) = cdn(3) - node_%z

                if(norm2(position(:)) < distance_min) then
                    distance_min = norm2(position(:))
                    rnodeID = nodeID
                    update = .true.
                end if

            end do check

        end do

    end function nearer_node

    logical function nearNode_check(cdn, nodeID, cubeID)
        real, intent(in) :: cdn(3)
        integer, intent(in) :: nodeID(3)
        integer, intent(in) :: cubeID
        real distance, position(3), edge(3)
        type(cube_inP3D) cube_
        type(node_inCUBE) node_

        cube_ = mesh(cubeID)

        node_ = cube_%nodes(nodeID(1), nodeID(2), nodeID(3))
        position(1) = cdn(1) - node_%x
        position(2) = cdn(2) - node_%y
        position(3) = cdn(3) - node_%z
        distance = norm2(position(:))

        edge(1) = cube_%nodes(2,1,1)%x - cube_%nodes(1,1,1)%x
        edge(2) = cube_%nodes(2,1,1)%y - cube_%nodes(1,1,1)%y
        edge(3) = cube_%nodes(2,1,1)%z - cube_%nodes(1,1,1)%z

        if(distance > 2.0*norm2(edge(:))) then
            nearNode_check = .false.
        else
            nearNode_check = .true.
        end if

    end function nearNode_check

    function get_velocity_f(nodeID, cubeID) result(velocity)
        integer, intent(in) :: nodeID(3)
        integer, intent(in) :: cubeID
        type(cube_inP3D) cube_
        real velocity(3)

        cube_ = mesh(cubeID)

        velocity(:) = cube_%nodes(nodeID(1),nodeID(2),nodeID(3))%f(2:4) &
                        / cube_%nodes(nodeID(1),nodeID(2),nodeID(3))%f(1)

    end function get_velocity_f

    subroutine deallocation_CUBE
        deallocate(mesh)
        deallocate(mesh_full)
    end subroutine deallocation_CUBE

end module CUBE_mod