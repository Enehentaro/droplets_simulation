program sortMain_vtk_ver
    use kdTree_m
    use path_operator_m
    use unstructuredGrid_mod
    implicit none 
    type(UnstructuredGrid) grid
    real, allocatable :: xyz(:,:)
    type(kdTree) kd_tree
    real droplet_position(3)
    integer nearest_ID
    integer n_unit 
    integer i, iimx, kkmx
    character(:), allocatable :: vtkFName
    character(10), parameter :: output_dir = 'test_check'
    real kdTree_startTime, kdTree_endTime, fullSearch_startTime, fullSearch_endTime
    integer nearestID, j
    real, allocatable :: distance(:)

    vtkFName = "sample.vtk"
            
    call grid%setupWithFlowFieldFile(vtkFName)

    iimx = size(grid%CELLs)
    kkmx = size(grid%NODEs)
    allocate(xyz(3, size(grid%CElls)))

    call make_directory(output_dir)

    do i = 1, iimx
        xyz(:, i) = grid%CELLs(i)%center(:)
    end do

    open(newunit = n_unit, file = output_dir//"/before.txt", status = 'replace')
        do i = 1, iimx
            ! write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(3(f12.5))') xyz(:, i)
        end do
    close(n_unit)

    kd_tree = kdTree_(xyz)
    call kd_tree%saveAsDOT(xyz, output_dir//'/kdTree.dot')

    call cpu_time(kdTree_startTime)
    do i = 1, iimx
        print*, "cellID =",i
        droplet_position(:) = xyz(:, i)
        call kd_tree%search(xyz, droplet_position, nearest_ID)
    end do
    call cpu_time(kdTree_endTime)

    call cpu_time(fullSearch_startTime)
        allocate(distance(iimx))
        do j = 1, iimx
            droplet_position(:) = xyz(:, j)
            do i = 1, iimx
                distance(i) = norm2(xyz(:, i) - droplet_position(:))
            end do
            nearestID = minloc(distance, dim = 1)
            print*, "nearestID =", nearestID 
        end do
    call cpu_time(fullSearch_endTime)

    print*, "kdTree_elapsedTime =",kdTree_endTime - kdTree_startTime, "sec"
    print*, "fullSearch_elapsedTime =",fullSearch_endTime - fullSearch_startTime, "sec"

end program sortMain_vtk_ver