module plot3d_operator
    implicit none

    type :: node_inCUBE
        real x, y, z
        real, allocatable :: f(:)
    end type

    type :: cube_inP3D
        type(node_inCUBE), allocatable :: nodes(:,:,:)
    end type
    
    contains

    function read_plot3d_multigrid(fName) result(mesh)
        character(*), intent(in) :: fName
        type(cube_inP3D), allocatable :: mesh(:)
        real :: min_max(6)
        integer n_unit, i,j,k, i_cube, num_cube
        integer, allocatable :: ni(:), nj(:), nk(:)

        print*, 'READ_multigrid:', fName

        open(newunit=n_unit , form='unformatted', file=fName, status='old')

            read(n_unit) num_cube   ;print*, 'num_cube=', num_cube
            allocate(mesh(num_cube))
            allocate(ni(num_cube), nj(num_cube), nk(num_cube))

            read(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), i_cube=1,num_cube)

            print*, 'ni=', minval(ni(:)), maxval(ni(:))
            print*, 'nj=', minval(nj(:)), maxval(ni(:))
            print*, 'nk=', minval(nk(:)), maxval(ni(:))

            do i_cube = 1, num_cube
                allocate(mesh(i_cube)%nodes(ni(i_cube), nj(i_cube), nk(i_cube)))
            end do

            do i_cube = 1, num_cube
                read(n_unit) &
                    (((mesh(i_cube)%nodes(i,j,k)%x, i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
                    (((mesh(i_cube)%nodes(i,j,k)%y, i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
                    (((mesh(i_cube)%nodes(i,j,k)%z, i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )

                    ! (((iblank(i,j,k,i_cube), i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )
            end do

        close(n_unit)

        min_max = minMax_coordinates(mesh)

        print*, 'min_cdn =', min_max(1:3)
        print*, 'max_cdn =', min_max(4:6)

    end function read_plot3d_multigrid

    function minMax_coordinates(mesh) result(min_max)
        type(cube_inP3D), intent(in) :: mesh(:)
        real min_max(6)
        integer i_cube, num_cube

        num_cube = size(mesh)

        min_max(1:3) = 1.e9
        min_max(4:6) = -1.e9

        do i_cube = 1, num_cube
            min_max(1) = min(min_max(1), minval(mesh(i_cube)%nodes%x))
            min_max(2) = min(min_max(2), minval(mesh(i_cube)%nodes%y))
            min_max(3) = min(min_max(3), minval(mesh(i_cube)%nodes%z))
            min_max(4) = max(min_max(4), maxval(mesh(i_cube)%nodes%x))
            min_max(5) = max(min_max(5), maxval(mesh(i_cube)%nodes%y))
            min_max(6) = max(min_max(6), maxval(mesh(i_cube)%nodes%z))
        end do

    end function minMax_coordinates

    subroutine write_plot3d_multigrid(mesh_in, fName)
        type(cube_inP3D), intent(in) :: mesh_in(:)
        character(*), intent(in) :: fName
        integer n_unit, i,j,k, i_cube, num_cube
        integer, allocatable :: ni(:), nj(:), nk(:)

        num_cube = size(mesh_in)

        allocate(ni(num_cube), nj(num_cube), nk(num_cube))

        do i_cube = 1, num_cube
            ni(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
            nj(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
            nk(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
        end do

        open(newunit=n_unit , form='unformatted', file=fName, status='replace')

            write(n_unit) num_cube

            write(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), i_cube=1,num_cube)

            do i_cube = 1, num_cube
                write(n_unit) &
                    ((( mesh_in(i_cube)%nodes(i,j,k)%x, &
                        i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
                    ((( mesh_in(i_cube)%nodes(i,j,k)%y, &
                        i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
                    ((( mesh_in(i_cube)%nodes(i,j,k)%z, &
                        i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )

                    ! (((iblank(i,j,k,i_cube), i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) )
            end do

        close(n_unit)
    end subroutine write_plot3d_multigrid

    subroutine write_plot3d_f(mesh_in, fName)
        type(cube_inP3D), intent(in) :: mesh_in(:)
        character(*), intent(in) :: fName
        integer n_unit, i,j,k,l, i_cube, num_cube
        integer, allocatable :: ni(:), nj(:), nk(:), nf(:)

        num_cube = size(mesh_in)

        allocate(ni(num_cube), nj(num_cube), nk(num_cube), nf(num_cube))

        do i_cube = 1, num_cube
            ni(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
            nj(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
            nk(i_cube) = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
        end do

        do i_cube = 1, num_cube
            nf(i_cube) = size(mesh_in(i_cube)%nodes(1,1,1)%f(:), dim=1)
        end do

        open(newunit=n_unit , form='unformatted', file=fName, status='replace')

            write(n_unit) num_cube

            write(n_unit) (ni(i_cube), nj(i_cube), nk(i_cube), nf(i_cube), i_cube=1,num_cube)

            do i_cube = 1, num_cube
                write(n_unit) &
                (((( mesh_in(i_cube)%nodes(i,j,k)%f(l), &
                i = 1, ni(i_cube) ), j = 1, nj(i_cube) ), k = 1, nk(i_cube) ), &
                l = 1, nf(i_cube))
            end do

        close(n_unit)
    end subroutine write_plot3d_f

    subroutine write_plot3d_asVTK(mesh_in, fName, n_cell)   !ひとつのcubeをひとつの節点とみなしてVTK出力
        type(cube_inP3D), intent(in) :: mesh_in(:)
        character(*), intent(in) :: fName
        integer, intent(in) :: n_cell    !CUBE一辺当たりのセル数
        integer n_unit, i_cube, num_cube, i,j,k, num_func, i_node, num_nodes, num_cells, i_cell
        integer i_max, j_max, k_max, delta_i, delta_j, delta_k, i1,i2,i3,i4, i_n
        real :: cdn(3), velocity(3)

        num_cube = size(mesh_in)
        num_nodes = num_cube*((n_cell + 1)**3)
        num_cells = num_cube*(n_cell**3)

        print*, 'VTK informations:'
        print*, 'nodes =', num_nodes
        print*, 'cells=', num_cells

        open(newunit=n_unit , form='formatted', file=fName, status='replace')

            write(n_unit, '(A)') '# vtk DataFile Version 2.0'
            write(n_unit, '(A)') 'Header'
            write(n_unit, '(A)') 'ASCII'
            write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

            write(n_unit, *) 'POINTS', num_nodes, 'float'
            do i_cube = 1, num_cube
                i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
                j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
                k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
                delta_i = (i_max - 1) / n_cell
                delta_j = (j_max - 1) / n_cell
                delta_k = (k_max - 1) / n_cell
                do k = 1, k_max, delta_k
                    do j = 1, j_max, delta_j
                        do i = 1, i_max, delta_i
                            cdn(1) = mesh_in(i_cube)%nodes(i,j,k)%x
                            cdn(2) = mesh_in(i_cube)%nodes(i,j,k)%y
                            cdn(3) = mesh_in(i_cube)%nodes(i,j,k)%z
                            write(n_unit, '(3(E15.8e2,2X))') cdn(:)
                        end do
                    end do
                end do
            end do

            write(n_unit, *) 'CELLS ', num_cells, num_cells*9
            do i_cube = 1, num_cube
                i_node =  (i_cube - 1)*((n_cell + 1)**3)    !CUBEの基準点ID
                i_n = 0
                do i_cell = 1, n_cell**3
                    i1 = i_node + i_n    !セルの基準点ID
                    i2 = i1 + n_cell + 1
                    i3 = i1 + (n_cell + 1)**2
                    i4 = i3 + n_cell + 1
                    write(n_unit, *)  8, i1, i1+1, i2, i2+1, i3, i3+1, i4, i4+1
                    i_n = i_n + 1
                    if(mod(i_n +1, n_cell +1) == 0) i_n = i_n + 1
                    if(mod(i_n +1 +n_cell, (n_cell +1)**2)==0) i_n = i_n + (n_cell + 1)
                end do
            end do
            write(n_unit, *) 'CELL_TYPES', num_cells
            do i_cell = 1, num_cells
                write(n_unit, *) 11 !11:直方体
            end do

            num_func = size(mesh_in(1)%nodes(1,1,1)%f(:))
            write(n_unit, *) 'POINT_DATA', num_nodes

            write(n_unit, '(A)') 'SCALARS density float'
            write(n_unit, '(A)') 'LOOKUP_TABLE default'
            do i_cube = 1, num_cube
                i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
                j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
                k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
                delta_i = (i_max - 1) / n_cell
                delta_j = (j_max - 1) / n_cell
                delta_k = (k_max - 1) / n_cell
                do k = 1, k_max, delta_k
                    do j = 1, j_max, delta_j
                        do i = 1, i_max, delta_i
                            write(n_unit, '(E15.8e2)') mesh_in(i_cube)%nodes(i,j,k)%f(1)
                        end do
                    end do
                end do
            end do

            write(n_unit, '(A)') 'VECTORS velocity float'
            do i_cube = 1, num_cube
                i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
                j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
                k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
                delta_i = (i_max - 1) / n_cell
                delta_j = (j_max - 1) / n_cell
                delta_k = (k_max - 1) / n_cell
                do k = 1, k_max, delta_k
                    do j = 1, j_max, delta_j
                        do i = 1, i_max, delta_i
                            velocity(:) = mesh_in(i_cube)%nodes(i,j,k)%f(2:4) / mesh_in(i_cube)%nodes(i,j,k)%f(1)
                            write(n_unit, '(3(E15.8e2,2X))') velocity(:)                                
                        end do
                    end do
                end do
            end do

        close(n_unit)
    end subroutine write_plot3d_asVTK

    subroutine read_plot3d_f(mesh, fName)
        type(cube_inP3D), intent(inout) :: mesh(:)
        character(*), intent(in) :: fName
        integer num_cube
        integer n_unit, i,j,k,l, i_cube
        real, allocatable :: min_f(:), max_f(:)
        integer, allocatable :: ni(:), nj(:), nk(:), nf(:)

        print*, 'READ_function:', fName

        open(newunit=n_unit , form='unformatted', file=fName, status='old')

            read(n_unit) num_cube   ;print*, 'num_cube=', num_cube
            if(num_cube /= size(mesh)) then
                print*, 'ERROR:the number of cube is not macth', num_cube, size(mesh)
                stop
            end if

            allocate(ni(num_cube), nj(num_cube), nk(num_cube), nf(num_cube))

            read(n_unit) (ni(i_cube),nj(i_cube),nk(i_cube),nf(i_cube), i_cube=1,num_cube)

            print*, 'nf=', minval(nf(:)), maxval(nf(:))

            if(.not. allocated(mesh(1)%nodes(1,1,1)%f)) then

                do i_cube = 1, num_cube
                    do k = 1, nk(i_cube)
                        do j = 1, nj(i_cube)
                            do i = 1, ni(i_cube)
                                allocate(mesh(i_cube)%nodes(i,j,k)%f(nf(i_cube)))
                            end do
                        end do
                    end do
                end do

            end if

            do i_cube = 1, num_cube
                read(n_unit) &
                    ((((mesh(i_cube)%nodes(i,j,k)%f(l), &
                        i = 1, ni(i_cube)), j = 1, nj(i_cube)), k = 1, nk(i_cube)), l = 1, nf(i_cube))
            end do

            allocate(min_f(maxval(nf(:))), source=1.e9)
            allocate(max_f(maxval(nf(:))), source=-1.e9)

            do i_cube = 1, num_cube
                do l = 1, nf(i_cube)
                    do k = 1, nk(i_cube)
                        do j = 1, nj(i_cube)
                            do i = 1, ni(i_cube)
                                min_f(l) = min(min_f(l), mesh(i_cube)%nodes(i,j,k)%f(l))
                                max_f(l) = max(max_f(l), mesh(i_cube)%nodes(i,j,k)%f(l))
                            end do
                        end do
                    end do
                end do
            end do

            print*, min_f
            print*, max_f

        close(n_unit)
    end subroutine read_plot3d_f

    function extract_cube(mesh_in, min_cdn, max_cdn) result(mesh_extracted)
        type(cube_inP3D), intent(in) :: mesh_in(:)
        type(cube_inP3D), allocatable :: mesh_extracted(:)
        real, intent(in) :: min_cdn(3), max_cdn(3)
        real :: center(3)
        integer i_cube, l, cube_cnt, num_cube
        logical extract
        integer :: original_ID(size(mesh_in))

        num_cube = size(mesh_in)

        cube_cnt = 0
         
        do i_cube = 1, num_cube
            extract = .false.

            center(:) = get_center_position(mesh_in(i_cube)%nodes)

            do l = 1, 3
                if(center(l) < min_cdn(l)) extract = .true.
                if(center(l) > max_cdn(l)) extract = .true.
            end do

            if(extract) cycle

            cube_cnt = cube_cnt +1
            original_ID(cube_cnt) = i_cube

        end do

        print*, 'cube_count=', cube_cnt

        allocate(mesh_extracted(cube_cnt))

        do i_cube = 1, cube_cnt
            mesh_extracted(i_cube) = mesh_in(original_ID(i_cube))
        end do

    end function extract_cube

    function change_resolution(mesh_in, resolution) result(mesh_changed)
        type(cube_inP3D), intent(in) :: mesh_in(:)
        type(cube_inP3D) :: mesh_changed(size(mesh_in))
        integer, intent(in) :: resolution   !一辺当たりのセル数
        integer i_cube, i,j,k, num_cube, i_max, j_max, k_max, delta_i, delta_j, delta_k
        integer ii,jj,kk

        num_cube = size(mesh_in)

        do i_cube = 1, num_cube
            allocate(mesh_changed(i_cube)%nodes(resolution+1, resolution+1, resolution+1))
            ! do k = 1, resolution+1
            !     do j = 1, resolution+1
            !         do i = 1, resolution+1
            !             nf = size(mesh_in(i_cube)%nodes(i,j,k)%f(:))
            !             allocate(mesh_changed(i_cube)%nodes(i,j,k)%f(nf))
            !         end do
            !     end do
            ! end do
        end do

        do i_cube = 1, num_cube
            i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
            j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
            k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
            delta_i = (i_max - 1) / resolution
            delta_j = (j_max - 1) / resolution
            delta_k = (k_max - 1) / resolution
            kk = 1
            do k = 1, k_max, delta_k
                jj = 1
                do j = 1, j_max, delta_j
                    ii = 1
                    do i = 1, i_max, delta_i
                        mesh_changed(i_cube)%nodes(ii,jj,kk) = mesh_in(i_cube)%nodes(i,j,k)
                        ii = ii + 1
                    end do
                    jj = jj + 1
                end do
                kk = kk + 1
            end do
        end do

    end function change_resolution

    function extract_function(mesh_in, nf) result(mesh_extracted)
        type(cube_inP3D), intent(in) :: mesh_in(:)
        type(cube_inP3D) :: mesh_extracted(size(mesh_in))
        integer, intent(in) :: nf   !関数の数
        integer i_cube, i,j,k, num_cube, i_max, j_max, k_max

        num_cube = size(mesh_in)

        do i_cube = 1, num_cube
            i_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=1)
            j_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=2)
            k_max = size(mesh_in(i_cube)%nodes(:,:,:), dim=3)
            allocate(mesh_extracted(i_cube)%nodes(i_max, j_max, k_max))
            do k = 1, k_max
                do j = 1, j_max
                    do i = 1, i_max
                        allocate(mesh_extracted(i_cube)%nodes(i,j,k)%f(nf))
                        mesh_extracted(i_cube)%nodes(i,j,k)%f(:) = mesh_in(i_cube)%nodes(i,j,k)%f(:nf)
                    end do
                end do
            end do
        end do

    end function extract_function

    function get_center_position(nodes) result(center)
        type(node_inCUBE) :: nodes(:,:,:)
        real :: center(3)
        integer ni,nj,nk,num_nodes

        ni = size(nodes(:,:,:), dim=1)
        nj = size(nodes(:,:,:), dim=2)
        nk = size(nodes(:,:,:), dim=3)

        num_nodes = ni * nj * nk
        center(1) = sum(nodes(:,:,:)%x)
        center(2) = sum(nodes(:,:,:)%y) 
        center(3) = sum(nodes(:,:,:)%z)
        center(:) = center(:) / num_nodes

    end function get_center_position

    real function mean_value(nodes, n_func)
        type(node_inCUBE) :: nodes(:,:,:)
        integer, intent(in) :: n_func
        integer ni,nj,nk,num_nodes, i,j,k

        ni = size(nodes(:,:,:), dim=1)
        nj = size(nodes(:,:,:), dim=2)
        nk = size(nodes(:,:,:), dim=3)

        num_nodes = ni * nj * nk

        mean_value = 0.0
        do k = 1, nk
            do j = 1, nj
                do i = 1, ni
                    mean_value = mean_value+ nodes(i,j,k)%f(n_func)
                end do
            end do
        end do

        mean_value = mean_value / num_nodes

    end function mean_value

end module plot3d_operator