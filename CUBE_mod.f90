module CUBE_mod
    use plot3d_operator
    use stl_reader
    implicit none

    type(cube_inP3D), allocatable, private :: mesh_full(:)
    type(cube_inP3D), allocatable, private :: mesh(:)

    type, extends(face_inSTL) :: face_t
        real AB(3), BC(3), CA(3)
    end type face_t

    type(face_t), allocatable :: faceShape(:)

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

    subroutine read_faceShape(path)
        implicit none
        type(face_inSTL), allocatable :: solid1(:), solid2(:), solid(:)
        character(*), intent(in) :: path
        character(len=99), allocatable :: stl_fname(:)
        integer i, num_face, num_file
    
        stl_fname = read_stl_fName()
        num_file = size(stl_fname)
        print*, 'num_stlFile=', num_file
        solid1 = read_stl(trim(path)//trim(stl_fname(1)))
        solid = solid1
        i = 1
        do while(i < num_file)
            if(allocated(solid)) deallocate(solid)

            i = i + 1
            solid2 = read_stl(trim(path)//trim(stl_fname(i)))

            solid = [solid1, solid2]

            deallocate(solid1, solid2)
            solid1 = solid

        end do

        num_face = size(solid)
        allocate(faceShape(num_face))
        faceShape%face_inSTL = solid
        print*, 'Total Faces on STL:', num_face

    end subroutine read_faceShape

    function read_stl_fName() result(stl_fnames)
        character(len=99), allocatable :: stl_fnames(:)
        integer n_unit, num_rec, ios, i
        character(len=99) str

        num_rec = 0

        open(newunit=n_unit, file='stl_fnames.txt', status='old')
            do        !レコード数を調べるループ
                read(n_unit, *, iostat=ios) str !ファイル終端であればiosに-1が返る
                if((trim(str) == '').or.(ios/=0)) exit    !空白か終端であればループ脱出
                num_rec = num_rec + 1    !カウント
            end do

            allocate(stl_fnames(num_rec))

            rewind(n_unit)

            do i = 1, num_rec
                read(n_unit, *) stl_fnames(i) 
            end do

        close(n_unit)

    end function read_stl_fName

    subroutine set_faceShape
        implicit none
        integer i, num_face
        type(face_inSTL) :: face_

        num_face = size(faceShape)

        do i = 1, num_face
            face_ = faceShape(i)%face_inSTL
            faceShape(i)%AB(:) = face_%nodes(2)%coordinate(:) - face_%nodes(1)%coordinate(:)
            faceShape(i)%BC(:) = face_%nodes(3)%coordinate(:) - face_%nodes(2)%coordinate(:)
            faceShape(i)%CA(:) = face_%nodes(1)%coordinate(:) - face_%nodes(3)%coordinate(:)
        end do

    end subroutine set_faceShape

    logical function adhesion_check_inSTL(X)
        integer i, num_face
        real, intent(in) :: X(3)
        real r_vec(3), inner, AS(3), BS(3), CS(3), Across(3), Bcross(3), Ccross(3)

        adhesion_check_inSTL = .false.

        num_face = size(faceShape)
        !$OMP parallel do private(r_vec, inner, AS, BS, CS, Across, Bcross, Ccross)
        do i = 1, num_face

            r_vec(:) = X(:) - faceShape(i)%nodes(1)%coordinate(:)  !点Aから飛沫への位置ベクトル
            inner = inner_product(r_vec, faceShape(i)%n_vector(:))    !位置ベクトルと法線ベクトルとの内積
            AS(:) = r_vec(:) - inner * faceShape(i)%n_vector(:)  !位置ベクトルを面へ投影
            BS(:) = - faceShape(i)%AB(:) + AS(:)
            CS(:) = faceShape(i)%CA(:) + AS(:)
            Across = cross_product(faceShape(i)%AB, AS)
            Bcross = cross_product(faceShape(i)%BC, BS)
            Ccross = cross_product(faceShape(i)%CA, CS)

            if(inner_product(Across, Bcross) < 0.0) cycle !三角形の外部
            if(inner_product(Across, Ccross) < 0.0) cycle

            if (abs(inner) <= 1.0d-2) then
                print*, 'inner=', inner
                adhesion_check_inSTL = .true.
            end if
        end do
        !$OMP end parallel do

    end function adhesion_check_inSTL

    function cross_product(a, b) result(cross)
        real,intent(in) :: a(3), b(3)
        real cross(3)

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function

    function inner_product(a, b) result(inner)
        real,intent(in) :: a(3), b(3)
        real inner

        inner = sum(a(:)*b(:))

    end function

end module CUBE_mod