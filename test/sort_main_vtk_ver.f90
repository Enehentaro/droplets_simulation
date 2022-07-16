program sortMain_vtk_ver
    use kdTree_m
    use path_operator_m
    use unstructuredGrid_m
    implicit none 
    type(UnstructuredGrid) grid
    real, allocatable :: xyz(:,:)
    type(kdTree) kd_tree
    real droplet_position(3)
    integer nearest_ID
    integer n_unit 
    integer i, iimx!, kkmx
    character(:), allocatable :: vtkFName
    character(10), parameter :: output_dir = 'test_check'

    vtkFName = "sample.vtk"
            
    call grid%setupWithFlowFieldFile(vtkFName)

    iimx = grid%get_info('cell')

    call make_directory(output_dir)

    xyz = grid%get_cellCenters()

    open(newunit = n_unit, file = output_dir//"/before.txt", status = 'replace')
        do i = 1, iimx
            ! write(n_unit,'(I3)', advance='no') before(i)%originID
            write(n_unit,'(3(f12.5))') xyz(:, i)
        end do
    close(n_unit)

    kd_tree = kdTree_(xyz)
    call kd_tree%saveAsDOT(xyz, output_dir//'/kdTree.dot')

    do i = 1, iimx
        print*, "cellID =",i
        droplet_position(:) = xyz(:, i)
        call kd_tree%search(xyz, droplet_position, nearest_ID)
    end do

end program sortMain_vtk_ver