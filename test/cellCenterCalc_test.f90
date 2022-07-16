!セル重心計算サブルーチンが2箇所に存在し、統一しようとしたものの絡まってしまったので断念
!ふたつのサブルーチン結果が一致するかをテストすることで齟齬が起きないよう対処（その場しのぎ）
program cellCenterCalc_test
    use unstructuredGrid_m
    use vtkMesh_operator_m
    use array_m
    implicit none
    type(UnstructuredGrid) grid
    type(vtkMesh) vtk_mesh
    ! character(17), parameter :: cellCenterFName = 'CellCenters.array'
    real, allocatable :: centers1(:,:), centers2(:,:), vtk_vel(:,:)
    integer i, imax

    call grid%setupWithFlowFieldFile('SAX/sax_flow.vtk')
    centers1 = grid%get_cellCenters()

    call vtk_mesh%read('SAX/sax_flow.vtk', cellVector=vtk_vel)
    centers2 = vtk_mesh%get_cellCenters()

    if(.not.all(centers1 == centers2)) error stop

    imax = grid%get_info('cell')
    do i = 1, imax
        ! print*, grid%get_flowVelocityInCELL(i), vtk_vel(:,i)
        if(.not.all(grid%get_flowVelocityInCELL(i) == vtk_vel(:,i))) error stop
    end do

    ! call output_cellCenter('SAX/')

    contains

    ! subroutine output_cellCenter(dir)
    !     character(*), intent(in) :: dir
    !     real, allocatable :: centers(:,:)
    !     integer num_cell

    !     num_cell = grid%get_info('cell')

    !     centers = grid%get_cellCenters()

    !     call output_2dArray_asBinary(dir//cellCenterFName, centers)
        
    ! end subroutine
    
end program cellCenterCalc_test