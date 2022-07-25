program kdTree_test
    !!kdTreeによる探索結果と、厳密なnearest探索結果が一致するかどうかをテスト
    use kdTree_m
    use path_operator_m
    use unstructuredGrid_m
    implicit none 
    type(FlowFieldUnstructuredGrid) grid
    real, allocatable :: xyz(:,:)
    type(kdTree) kd_tree
    ! integer n_unit 
    integer iimx
    character(:), allocatable :: vtkFName, kd_treeFName
    character(10), parameter :: output_dir = 'test_check'
    logical existance

    vtkFName = "sample.vtk"
    kd_treeFName = "kdTreeOfsample.txt"
    call make_directory(output_dir)
    call grid%setupWithFlowFieldFile(vtkFName)
    iimx = grid%get_info('cell')
    xyz = grid%get_allOfCellCenters()



    inquire(file = output_dir//'/'//kd_treeFName, exist=existance)
    if(.not.existance) then
                
        !これいる？
        ! open(newunit = n_unit, file = output_dir//"/before.txt", status = 'replace')
        !     do i = 1, iimx
        !         ! write(n_unit,'(I3)', advance='no') before(i)%originID
        !         write(n_unit,'(3(f12.5))') xyz(:, i)
        !     end do
        ! close(n_unit)

        kd_tree = kdTree_(xyz)
        call kd_tree%saveAsDOT(xyz, output_dir//'/kdTree.dot')
        call kd_tree%saveAsTXT(output_dir//'/'//kd_treeFName)
    
    else

        call kd_tree%read_kdTree(output_dir//'/'//kd_treeFName, iimx)

    end if

    call test

    contains

    subroutine test
        !!乱数で発生させた点に対して、kdTreeによる探索結果と、厳密なnearest探索結果が一致するかどうかをテスト
        integer, parameter :: num_test = 10000
        integer n
        integer result_exact, result_kdTree
        real, dimension(3) :: min_cdn, max_cdn, delta, rand,  point

        call grid%get_MinMaxOfGrid(min_cdn, max_cdn)
        delta = max_cdn - min_cdn

        do n = 1, num_test
            call random_number(rand)
            point = min_cdn + delta*rand
            ! print*, point

            result_exact = grid%exact_nearest_search(point)
            call kd_tree%search(xyz, point, result_kdTree)

            if(result_exact /= result_kdTree) then
                print*, point
                print*, result_exact, xyz(:, result_exact)
                print*, result_kdTree, xyz(:, result_kdTree)
                error stop
            end if

        end do

    end subroutine

end program kdTree_test