program sort_main_vtk_ver
    use kdTree_m
    use path_operator_m
    use unstructuredGrid_m
    implicit none 
    type(FlowFieldUnstructuredGrid) grid
    real, allocatable :: xyz(:,:)
    type(kdTree) kd_tree
    real droplet_position(3)
    integer nearest_ID
    integer n_unit 
    integer i, iimx
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
                
        open(newunit = n_unit, file = output_dir//"/before.txt", status = 'replace')
            do i = 1, iimx
                ! write(n_unit,'(I3)', advance='no') before(i)%originID
                write(n_unit,'(3(f12.5))') xyz(:, i)
            end do
        close(n_unit)

        kd_tree = kdTree_(xyz)
        call kd_tree%saveAsDOT(xyz, output_dir//'/kdTree.dot')
        call kd_tree%saveAsTXT(output_dir//'/'//kd_treeFName)
    
        do i = 1, iimx
            print*, "cellID =",i
            droplet_position(:) = xyz(:, i)
            call kd_tree%search(xyz, droplet_position, nearest_ID)
            print*, "nearest_ID =", nearest_ID
        end do

    else

        call kd_tree%read_kdTree(output_dir//'/'//kd_treeFName, iimx)
    
        do i = 1, iimx
            print*, "cellID =",i
            droplet_position(:) = xyz(:, i)
            call kd_tree%search(xyz, droplet_position, nearest_ID)
            print*, "nearest_ID =", nearest_ID
        end do

    end if

end program sort_main_vtk_ver