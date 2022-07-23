!>セル重心計算がうまく行っているかをテストする。
!>具体的には、重心をテトラの内外判定にかけている。
program cellCenterCalc_test
    use unstructuredGrid_m
    use geometry_m
    implicit none
    type(FlowFieldUnstructuredGrid) grid
    ! character(17), parameter :: cellCenterFName = 'CellCenters.array'
    real, allocatable :: centers(:,:), vertices(:,:)
    real center(3)
    integer i, imax
    integer n

    call grid%setupWithFlowFieldFile('SAX/sax_flow.vtk')    !このサブルーチン内で重心計算も行われる

    imax = grid%get_info('cell')
    centers = grid%get_allOfCellCenters()

    do i = 1, imax
        vertices = grid%get_cellVerticesOf(i)
        center = centers(:,i)

        !重心がテトラ内部になければエラー
        if(.not.insideJudgment_tetra(vertices, center)) then
            print'("============================")'
            block
                real vol_sum, volume
                call insideJudgment_tetra_check(vertices, center, vol_sum, volume)
                print*, vol_sum, volume
            end block
            print'("============================")'
            print*, i, imax
            print'("============================")'
            do n = 1, size(vertices, dim=2)
                print*, vertices(:,n)
            end do
            print'("============================")'
            print*, center
            print'("============================")'
            error stop
        end if

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