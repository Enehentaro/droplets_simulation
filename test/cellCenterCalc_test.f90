!>セル重心計算がうまく行っているかをテストする。
!>具体的には、重心をテトラの内外判定にかけている。
program cellCenterCalc_test
    use unstructuredGrid_m
    use geometry_m
    implicit none
    type(FlowFieldUnstructuredGrid) grid
    ! character(17), parameter :: cellCenterFName = 'CellCenters.array'
    real, allocatable :: centers(:,:), vertices(:,:)
    integer i, imax

    call grid%setupWithFlowFieldFile('SAX/sax_flow.vtk')    !このサブルーチン内で重心計算も行われる
    centers = grid%get_allOfCellCenters()

    do i = 1, imax
        vertices = grid%get_cellVerticesOf(i)
        !重心がテトラ内部になければエラー
        if(.not.insideJudgment_tetra(vertices, centers(:,i))) error stop
    end do

    ! contains

    ! subroutine output_cellCenter(dir)
    !     use array_m
    !     character(*), intent(in) :: dir
    !     real, allocatable :: centers(:,:)
    !     integer num_cell

    !     num_cell = grid%get_info('cell')

    !     centers = grid%get_cellCenters()

    !     call output_2dArray_asBinary(dir//cellCenterFName, centers)
        
    ! end subroutine
    
end program cellCenterCalc_test