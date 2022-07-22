!>CUBE出力ファイル形式：Plot3Dを取り扱うためのモジュール
module plot3d_operator
    implicit none
    private

    !>エリア構造体
    type area_t
        real, dimension(3) :: min, max  !最小座標と最大座標
    end type

    !>cube（Plot3D形式における立方体）構造体
    type cube_inP3D
        real, allocatable :: nodes(:,:,:,:)
            !!節点座標配列（i,j,k, xyz）

        real, allocatable :: f(:,:,:,:) 
            !!保存量配列（i,j,k, f）

        type(area_t) area
        contains
        procedure, public :: isIncluded
        procedure areaOfCube, nearest_nodeID!, nearer_nodeID
    end type

    !>節点情報構造体
    type, public :: plot3dNodeInfo  
        integer cubeID, nodeID(3)   !節点が属するcubeID、そのcube内の節点ID（i,j,kそれぞれのIDを指定するので要素数3）
    end type

    !>Plot3Dメッシュクラス
    type, public :: Plot3dMesh  
        private
        type(cube_inP3D), allocatable :: cubes(:)
            !!cube配列

        contains
        procedure areaOfMesh
        procedure, public :: read_plot3d_function, get_numCube, get_velocity, get_cubeShape
        procedure, public :: get_cubeID_contains, nearestNodeInfo
    end type

    public :: read_plot3d_multigrid !Plot3Dメッシュクラスの事実上コンストラクタ
    
    contains

    !>メッシュファイル（.g）を読み込んでメッシュクラスを返す関数
    function read_plot3d_multigrid(fName) result(mesh)
        character(*), intent(in) :: fName
        type(Plot3dMesh) mesh
        integer n_unit, i,j,k, i_cube, num_cube, xyz
        integer, allocatable :: ni(:), nj(:), nk(:)

        print*, 'READ_multigrid:', fName

        open(newunit=n_unit , form='unformatted', file=fName, status='old', action='read')

            read(n_unit) num_cube   ;print*, 'num_cube=', num_cube
            allocate(mesh%cubes(num_cube))
            allocate(ni(num_cube), nj(num_cube), nk(num_cube))

            read(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), i_cube=1,num_cube)

            print*, 'ni=', minval(ni(:)), maxval(ni(:))
            print*, 'nj=', minval(nj(:)), maxval(ni(:))
            print*, 'nk=', minval(nk(:)), maxval(ni(:))

            do i_cube = 1, num_cube
                allocate(mesh%cubes(i_cube)%nodes(ni(i_cube), nj(i_cube), nk(i_cube), 3))
            end do

            do i_cube = 1, num_cube
                read(n_unit) ((((   mesh%cubes(i_cube)%nodes(i,j,k, xyz), &
                    i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), xyz = 1,3)
                    ! (((iblank(i,j,k,i_cube), i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )

                mesh%cubes(i_cube)%area = mesh%cubes(i_cube)%areaOfCube()
            end do

        close(n_unit)

        block
            type(area_t) area

            area = mesh%areaOfMesh()
            print*, 'min_cdn =', area%min
            print*, 'max_cdn =', area%max
        end block

    end function

    !>メッシュのエリアを構造体で返す関数
    type(area_t) function areaOfMesh(self)
        class(Plot3dMesh), intent(in) :: self
        integer i_cube, num_cube

        num_cube = size(self%cubes)

        areaOfMesh%min(:) = 1.e9
        areaOfMesh%max(:) = -1.e9

        do i_cube = 1, num_cube
            areaOfMesh%min(1) = min(areaOfMesh%min(1), self%cubes(i_cube)%area%min(1))
            areaOfMesh%min(2) = min(areaOfMesh%min(2), self%cubes(i_cube)%area%min(2))
            areaOfMesh%min(3) = min(areaOfMesh%min(3), self%cubes(i_cube)%area%min(3))
            areaOfMesh%max(1) = max(areaOfMesh%max(1), self%cubes(i_cube)%area%max(1))
            areaOfMesh%max(2) = max(areaOfMesh%max(2), self%cubes(i_cube)%area%max(2))
            areaOfMesh%max(3) = max(areaOfMesh%max(3), self%cubes(i_cube)%area%max(3))
        end do

    end function

    !>cubeのエリアを構造体で返す関数
    type(area_t) function areaOfCube(self)
        class(cube_inP3D), intent(in) :: self

        areaOfCube%min(1) = minval(self%nodes(:,:,:,1))
        areaOfCube%min(2) = minval(self%nodes(:,:,:,2))
        areaOfCube%min(3) = minval(self%nodes(:,:,:,3))
        areaOfCube%max(1) = maxval(self%nodes(:,:,:,1))
        areaOfCube%max(2) = maxval(self%nodes(:,:,:,2))
        areaOfCube%max(3) = maxval(self%nodes(:,:,:,3))

    end function

    !>cubeのエリア内に任意座標（引数）が含まれているかを返す関数
    logical function isIncluded(self, cdn)
        class(cube_inP3D), intent(in) :: self
        real, intent(in) :: cdn(3)

        isIncluded = .false.

        if ((cdn(1) >= self%area%min(1)).and.(cdn(2) >= self%area%min(2)).and.(cdn(3) >= self%area%min(3))) then

            if  ((cdn(1) <= self%area%max(1)).and.(cdn(2) <= self%area%max(2)).and.(cdn(3) <= self%area%max(3)))then

                isIncluded = .true.

            end if

        end if

    end function

    !>任意座標（引数）を含むcubeをメッシュの中から探してそのIDを返す関数
    integer function get_cubeID_contains(self, cdn)
        class(Plot3dMesh), intent(in) :: self
        real, intent(in) :: cdn(3)
        integer i_cube, num_cube

        get_cubeID_contains = 0
        num_cube = self%get_numCube()

        do i_cube = 1, num_cube
            if(self%cubes(i_cube)%isIncluded(cdn)) then
                get_cubeID_contains = i_cube
                return
            end if
        end do

        print*, 'ERROR : The Point is Out of Area.', cdn
        error stop

    end function

    !>任意座標（引数）に最近傍な節点を探してそのIDを返す関数
    function nearest_nodeID(self, cdn) result(nodeID)
        class(cube_inP3D), intent(in) :: self
        real, intent(in) :: cdn(3)
        integer nodeID(3)
        integer i,j,k, cubeShape(4)
        real relation(3)
        real, allocatable :: distance(:,:,:)

        cubeShape = shape(self%nodes)
        allocate(distance(cubeShape(1), cubeShape(2), cubeShape(3)), source=1.e9)
        !$OMP parallel private(relation)
        !$OMP do collapse(3)
        do k = 1, cubeShape(3)
            do j = 1, cubeShape(2)
                do i = 1, cubeShape(1)
                    relation(:) = cdn(:) - self%nodes(i,j,k, :)
                    distance(i,j,k) = norm2(relation(:))
                end do
            end do
        end do
        !$OMP end do
        !$OMP end parallel

        nodeID = minloc(distance)

    end function

    ! function nearer_nodeID(self, cdn, inodeID) result(rnodeID)
    !     class(cube_inP3D), intent(in) :: self
    !     real, intent(in) :: cdn(3)
    !     integer, intent(in) :: inodeID(3)
    !     integer nodeID(3), rnodeID(3), i, maxID
    !     real distance_min, relation(3)
    !     integer,parameter :: diff(3,6) &
    !                         = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], shape(diff))
    !     logical update
    !     type(cube_inP3D) cube_
    !     type(node_inCUBE) node_

    !     node_ = self%nodes(inodeID(1), inodeID(2), inodeID(3))
    !     relation(1) = cdn(1) - node_%x
    !     relation(2) = cdn(2) - node_%y
    !     relation(3) = cdn(3) - node_%z
    !     distance_min = norm2(relation(:))
    !     rnodeID = inodeID
    !     update = .true.

    !     maxID = size(cube_%nodes(:,:,:), dim=1)

    !     do while(update)
    !         update = .false.

    !         check : do i = 1, 6

    !             nodeID(:) = rnodeID(:) + diff(:,i)

    !             if((minval(nodeID) < 1).or.(maxval(nodeID) > maxID)) cycle check

    !             node_ = cube_%nodes(nodeID(1), nodeID(2), nodeID(3))
    !             relation(1) = cdn(1) - node_%x
    !             relation(2) = cdn(2) - node_%y
    !             relation(3) = cdn(3) - node_%z

    !             if(norm2(relation(:)) < distance_min) then
    !                 distance_min = norm2(relation(:))
    !                 rnodeID = nodeID
    !                 update = .true.
    !             end if

    !         end do check

    !     end do

    ! end function

    !>任意座標（引数）に最近傍な節点を探し、その情報（cubeID,nodeID）を返す関数
    type(plot3dNodeInfo) function nearestNodeInfo(self, cdn)
        class(Plot3dMesh), intent(in) :: self
        real, intent(in) :: cdn(3)
        integer cubeID!, nodeID(3)

        cubeID = self%get_cubeID_contains(cdn)
        
        nearestNodeInfo%cubeID = cubeID
        nearestNodeInfo%nodeID(:) = self%cubes(cubeID)%nearest_nodeID(cdn)

        ! nodeID = self%cubes(cubeID)%nearest_nodeID(cdn, 4)
        ! nearestNodeInfo%nodeID(:) = self%cubes(cubeID)%nearer_nodeID(cdn, nodeID)

    end function

    ! subroutine write_plot3d_multigrid(mesh_in, fName)
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     character(*), intent(in) :: fName
    !     integer n_unit, i,j,k, i_cube, num_cube
    !     integer, allocatable :: ni(:), nj(:), nk(:)

    !     num_cube = size(mesh_in)

    !     allocate(ni(num_cube), nj(num_cube), nk(num_cube))

    !     do i_cube = 1, num_cube
    !         ni(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !         nj(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !         nk(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !     end do

    !     open(newunit=n_unit , form='unformatted', file=fName, status='replace')

    !         write(n_unit) num_cube

    !         write(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), i_cube=1,num_cube)

    !         do i_cube = 1, num_cube
    !             write(n_unit) &
    !                 ((( mesh_in(i_cube)%nodes(i,j,k)%x, &
    !                     i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
    !                 ((( mesh_in(i_cube)%nodes(i,j,k)%y, &
    !                     i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
    !                 ((( mesh_in(i_cube)%nodes(i,j,k)%z, &
    !                     i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )

    !                 ! (((iblank(i,j,k,i_cube), i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )
    !         end do

    !     close(n_unit)
    ! end subroutine

    ! subroutine write_plot3d_f(mesh_in, fName)
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     character(*), intent(in) :: fName
    !     integer n_unit, i,j,k,l, i_cube, num_cube
    !     integer, allocatable :: ni(:), nj(:), nk(:), nf(:)

    !     num_cube = size(mesh_in)

    !     allocate(ni(num_cube), nj(num_cube), nk(num_cube), nf(num_cube))

    !     do i_cube = 1, num_cube
    !         ni(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !         nj(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !         nk(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !     end do

    !     do i_cube = 1, num_cube
    !         nf(i_cube) = size(mesh_in(i_cube)%nodes(1,1,1)%f(:), dim=1)
    !     end do

    !     open(newunit=n_unit , form='unformatted', file=fName, status='replace')

    !         write(n_unit) num_cube

    !         write(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), nf(i_cube), i_cube=1,num_cube)

    !         do i_cube = 1, num_cube
    !             write(n_unit) &
    !             (((( mesh_in(i_cube)%nodes(i,j,k)%f(l), &
    !             i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
    !             l = 1, nf(i_cube))
    !         end do

    !     close(n_unit)
    ! end subroutine

    ! subroutine write_plot3d_asVTK(mesh_in, fName, n_cell)   !ひとつのcubeをひとつの節点とみなしてVTK出力
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     character(*), intent(in) :: fName
    !     integer, intent(in) :: n_cell    !CUBE一辺当たりのセル数
    !     integer n_unit, i_cube, num_cube, i,j,k, num_func, i_node, num_nodes, num_cells, i_cell
    !     integer i_max, j_max, k_max, delta_i, delta_j, delta_k, i1,i2,i3,i4, i_n
    !     real :: cdn(3), velocity(3)

    !     num_cube = size(mesh_in)
    !     num_nodes = num_cube*((n_cell + 1)**3)
    !     num_cells = num_cube*(n_cell**3)

    !     print*, 'VTK informations:'
    !     print*, 'nodes =', num_nodes
    !     print*, 'cells=', num_cells

    !     open(newunit=n_unit , form='formatted', file=fName, status='replace')

    !         write(n_unit, '(A)') '# vtk DataFile Version 2.0'
    !         write(n_unit, '(A)') 'Header'
    !         write(n_unit, '(A)') 'ASCII'
    !         write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

    !         write(n_unit, *) 'POINTS', num_nodes, 'float'
    !         do i_cube = 1, num_cube
    !             i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !             j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !             k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !             delta_i = (i_max - 1) / n_cell
    !             delta_j = (j_max - 1) / n_cell
    !             delta_k = (k_max - 1) / n_cell
    !             do k = 1, k_max, delta_k
    !                 do j = 1, j_max, delta_j
    !                     do i = 1, i_max, delta_i
    !                         cdn(1) = mesh_in(i_cube)%nodes(i,j,k)%x
    !                         cdn(2) = mesh_in(i_cube)%nodes(i,j,k)%y
    !                         cdn(3) = mesh_in(i_cube)%nodes(i,j,k)%z
    !                         write(n_unit, '(3(E15.8e2,2X))') cdn(:)
    !                     end do
    !                 end do
    !             end do
    !         end do

    !         write(n_unit, *) 'CELLS ', num_cells, num_cells*9
    !         do i_cube = 1, num_cube
    !             i_node =  (i_cube - 1)*((n_cell + 1)**3)    !CUBEの基準点ID
    !             i_n = 0
    !             do i_cell = 1, n_cell**3
    !                 i1 = i_node + i_n    !セルの基準点ID
    !                 i2 = i1 + n_cell + 1
    !                 i3 = i1 + (n_cell + 1)**2
    !                 i4 = i3 + n_cell + 1
    !                 write(n_unit, *)  8, i1, i1+1, i2, i2+1, i3, i3+1, i4, i4+1
    !                 i_n = i_n + 1
    !                 if(mod(i_n +1, n_cell +1) == 0) i_n = i_n + 1
    !                 if(mod(i_n +1 +n_cell, (n_cell +1)**2)==0) i_n = i_n + (n_cell + 1)
    !             end do
    !         end do
    !         write(n_unit, *) 'CELL_TYPES', num_cells
    !         do i_cell = 1, num_cells
    !             write(n_unit, *) 11 !11:直方体
    !         end do

    !         num_func = size(mesh_in(1)%nodes(1,1,1)%f(:))
    !         write(n_unit, *) 'POINT_DATA', num_nodes

    !         write(n_unit, '(A)') 'SCALARS density float'
    !         write(n_unit, '(A)') 'LOOKUP_TABLE default'
    !         do i_cube = 1, num_cube
    !             i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !             j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !             k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !             delta_i = (i_max - 1) / n_cell
    !             delta_j = (j_max - 1) / n_cell
    !             delta_k = (k_max - 1) / n_cell
    !             do k = 1, k_max, delta_k
    !                 do j = 1, j_max, delta_j
    !                     do i = 1, i_max, delta_i
    !                         write(n_unit, '(E15.8e2)') mesh_in(i_cube)%nodes(i,j,k)%f(1)
    !                     end do
    !                 end do
    !             end do
    !         end do

    !         write(n_unit, '(A)') 'VECTORS velocity float'
    !         do i_cube = 1, num_cube
    !             i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !             j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !             k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !             delta_i = (i_max - 1) / n_cell
    !             delta_j = (j_max - 1) / n_cell
    !             delta_k = (k_max - 1) / n_cell
    !             do k = 1, k_max, delta_k
    !                 do j = 1, j_max, delta_j
    !                     do i = 1, i_max, delta_i
    !                         velocity(:) = mesh_in(i_cube)%nodes(i,j,k)%f(2:4) / mesh_in(i_cube)%nodes(i,j,k)%f(1)
    !                         write(n_unit, '(3(E15.8e2,2X))') velocity(:)                                
    !                     end do
    !                 end do
    !             end do
    !         end do

    !     close(n_unit)
    ! end subroutine

    !>メッシュクラスメソッド
    !>保存量ファイル（.f）を読み込み、メッシュクラスに格納する
    subroutine read_plot3d_function(self, fName)
        class(Plot3dMesh) self
        character(*), intent(in) :: fName
        integer num_cube
        integer n_unit, i,j,k,l, i_cube
        real, allocatable :: min_f(:), max_f(:)
        integer, allocatable :: ni(:), nj(:), nk(:), nf(:)

        print*, 'READ_function:', fName

        open(newunit=n_unit , form='unformatted', file=fName, status='old', action='read')

            read(n_unit) num_cube   ;print*, 'num_cube=', num_cube
            if(num_cube /= size(self%cubes)) then
                print*, 'ERROR:the number of cube is not macth', num_cube, size(self%cubes)
                error stop
            end if

            allocate(ni(num_cube), nj(num_cube), nk(num_cube), nf(num_cube))

            read(n_unit) (ni(i_cube),nj(i_cube),nk(i_cube),nf(i_cube), i_cube=1,num_cube)

            print*, 'nf=', minval(nf(:)), maxval(nf(:))

            if(.not.allocated(self%cubes(1)%f)) then
                do i_cube = 1, num_cube
                    allocate(self%cubes(i_cube)%f(ni(i_cube), nj(i_cube), nk(i_cube), nf(i_cube)))
                end do
            end if

            do i_cube = 1, num_cube
                read(n_unit) &
                    ((((self%cubes(i_cube)%f(i,j,k,l), &
                        i = 1, ni(i_cube)), j = 1, nj(i_cube)), k = 1, nk(i_cube)), l = 1, nf(i_cube))
            end do

            allocate(min_f(maxval(nf(:))), source=1.e9)
            allocate(max_f(maxval(nf(:))), source=-1.e9)

            do i_cube = 1, num_cube
                do l = 1, nf(i_cube)
                    min_f(l) = min(min_f(l), minval(self%cubes(i_cube)%f(:,:,:,l)))
                    max_f(l) = max(max_f(l), maxval(self%cubes(i_cube)%f(:,:,:,l)))
                end do
            end do

            print*, min_f
            print*, max_f

        close(n_unit)

    end subroutine

    ! function extract_cube(mesh_in, min_cdn, max_cdn) result(mesh_extracted)
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     type(cube_inP3D), allocatable :: mesh_extracted(:)
    !     real, intent(in) :: min_cdn(3), max_cdn(3)
    !     real :: center(3)
    !     integer i_cube, l, cube_cnt, num_cube
    !     logical extract
    !     integer :: original_ID(size(mesh_in))

    !     num_cube = size(mesh_in)

    !     cube_cnt = 0
         
    !     do i_cube = 1, num_cube
    !         extract = .false.

    !         center(:) = get_center_position(mesh_in(i_cube)%nodes)

    !         do l = 1, 3
    !             if(center(l) < min_cdn(l)) extract = .true.
    !             if(center(l) > max_cdn(l)) extract = .true.
    !         end do

    !         if(extract) cycle

    !         cube_cnt = cube_cnt +1
    !         original_ID(cube_cnt) = i_cube

    !     end do

    !     print*, 'cube_count=', cube_cnt

    !     allocate(mesh_extracted(cube_cnt))

    !     do i_cube = 1, cube_cnt
    !         mesh_extracted(i_cube) = mesh_in(original_ID(i_cube))
    !     end do

    ! end function

    ! function change_resolution(mesh_in, resolution) result(mesh_changed)
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     type(cube_inP3D) :: mesh_changed(size(mesh_in))
    !     integer, intent(in) :: resolution   !一辺当たりのセル数
    !     integer i_cube, i,j,k, num_cube, i_max, j_max, k_max, delta_i, delta_j, delta_k
    !     integer ii,jj,kk

    !     num_cube = size(mesh_in)

    !     do i_cube = 1, num_cube
    !         allocate(mesh_changed(i_cube)%nodes(resolution+1, resolution+1, resolution+1))
    !         ! do k = 1, resolution+1
    !         !     do j = 1, resolution+1
    !         !         do i = 1, resolution+1
    !         !             nf = size(mesh_in(i_cube)%nodes(i,j,k)%f(:))
    !         !             allocate(mesh_changed(i_cube)%nodes(i,j,k)%f(nf))
    !         !         end do
    !         !     end do
    !         ! end do
    !     end do

    !     do i_cube = 1, num_cube
    !         i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !         j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !         k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !         delta_i = (i_max - 1) / resolution
    !         delta_j = (j_max - 1) / resolution
    !         delta_k = (k_max - 1) / resolution
    !         kk = 1
    !         do k = 1, k_max, delta_k
    !             jj = 1
    !             do j = 1, j_max, delta_j
    !                 ii = 1
    !                 do i = 1, i_max, delta_i
    !                     mesh_changed(i_cube)%nodes(ii,jj,kk) = mesh_in(i_cube)%nodes(i,j,k)
    !                     ii = ii + 1
    !                 end do
    !                 jj = jj + 1
    !             end do
    !             kk = kk + 1
    !         end do
    !     end do

    ! end function

    ! function extract_function(mesh_in, nf) result(mesh_extracted)
    !     type(cube_inP3D), intent(in) :: mesh_in(:)
    !     type(cube_inP3D) :: mesh_extracted(size(mesh_in))
    !     integer, intent(in) :: nf   !関数の数
    !     integer i_cube, i,j,k, num_cube, i_max, j_max, k_max

    !     num_cube = size(mesh_in)

    !     do i_cube = 1, num_cube
    !         i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
    !         j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
    !         k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
    !         allocate(mesh_extracted(i_cube)%nodes(i_max, j_max, k_max))
    !         do k = 1, k_max
    !             do j = 1, j_max
    !                 do i = 1, i_max
    !                     allocate(mesh_extracted(i_cube)%nodes(i,j,k)%f(nf))
    !                     mesh_extracted(i_cube)%nodes(i,j,k)%f(:) = mesh_in(i_cube)%nodes(i,j,k)%f(:nf)
    !                 end do
    !             end do
    !         end do
    !     end do

    ! end function

    ! function get_center_position(nodes) result(center)
    !     type(node_inCUBE) :: nodes(:,:,:)
    !     real :: center(3)
    !     integer ni,nj,nk,num_nodes

    !     ni = size(nodes(:,:,:), dim=1)
    !     nj = size(nodes(:,:,:), dim=2)
    !     nk = size(nodes(:,:,:), dim=3)

    !     num_nodes = ni * nj * nk
    !     center(1) = sum(nodes(:,:,:)%x)
    !     center(2) = sum(nodes(:,:,:)%y) 
    !     center(3) = sum(nodes(:,:,:)%z)
    !     center(:) = center(:) / num_nodes

    ! end function

    ! real function mean_value(nodes, n_func)
    !     type(node_inCUBE) :: nodes(:,:,:)
    !     integer, intent(in) :: n_func
    !     integer ni,nj,nk,num_nodes, i,j,k

    !     ni = size(nodes(:,:,:), dim=1)
    !     nj = size(nodes(:,:,:), dim=2)
    !     nk = size(nodes(:,:,:), dim=3)

    !     num_nodes = ni * nj * nk

    !     mean_value = 0.0
    !     do k = 1, nk
    !         do j = 1, nj
    !             do i = 1, ni
    !                 mean_value = mean_value+ nodes(i,j,k)%f(n_func)
    !             end do
    !         end do
    !     end do

    !     mean_value = mean_value / num_nodes

    ! end function

    !>メッシュを構成するcube数を返す関数
    integer function get_numCube(self)
        class(Plot3dMesh), intent(in) :: self

        get_numCube = size(self%cubes)

    end function

    !>cubeの形状（i,j,k節点数）を返す関数
    function get_cubeShape(self) result(cubeShape)
        class(Plot3dMesh), intent(in) :: self
        integer shapeArray(4), cubeShape(3)

        shapeArray = shape(self%cubes(1)%nodes)
        cubeShape = shapeArray(1:3)

    end function

    !>任意節点（引数）における流速を返す関数
    function get_velocity(self, node) result(velocity)
        class(Plot3dMesh) self
        type(plot3dNodeInfo), intent(in) :: node
            !!節点
        
        type(cube_inP3D) cube
        real velocity(3)

        cube = self%cubes(node%cubeID)

        !質量流束（第2〜第4成分）を密度（第1成分）で割ることで流速を取得する
        velocity(:) = cube%f(node%nodeID(1),node%nodeID(2),node%nodeID(3), 2:4) &
                        / cube%f(node%nodeID(1),node%nodeID(2),node%nodeID(3), 1)

    end function

end module plot3d_operator
