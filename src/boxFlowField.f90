module boxOperator_m
    use boxCounter_m
    implicit none

    private

    type(boxCounter), allocatable, public :: box_array(:)

    type, public :: boxResult_t
        real flowVelocity(3)
    end type

    public setup_boxArray, output_boxVTK, output_flowCSV

    contains

    subroutine setup_boxArray
        box_array = get_box_array('.', 0)
    end subroutine

    subroutine output_flowCSV(fname, bResult)
        character(*), intent(in) :: fname
        type(boxResult_t), intent(in) :: bResult(:)
        integer n_unit, i

        print*, 'output: ', fname

        open(newunit=n_unit, file=fname, status='replace')
        
            write(n_unit, '("x,y,z,u,v,w")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%flowVelocity
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK(fname, bResult)
        use vtkMesh_operator_m
        character(*), intent(in) :: fname
        type(boxResult_t), intent(in) :: bResult(:)
        type(vtkMesh) boxMesh
        integer i, j, k, num_box
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
        real velArray(3, size(box_array))
        real, allocatable :: xyz(:,:)
        integer, allocatable :: vertices(:,:), types(:)
                                        
        num_box = size(box_array)
        allocate(xyz(3, num_box*8))
        allocate(vertices(8, num_box), types(num_box))
        do i = 1, num_box

            do j = 1, 8
                k = j + 8*(i-1)
                xyz(:,k) = box_array(i)%min_cdn(:) + box_array(i)%width(:)*trans(:,j)
                vertices(j,i) = k
            end do
            
            types(i) = 11

        end do

        do i = 1, num_box

            velArray(:,i) = bResult(i)%flowVelocity

        end do

        boxMesh = vtkMesh_(xyz, vertices, types)

        call boxMesh%output(fname, cellVector=velArray, vectorName='VEL')

    end subroutine

end module boxOperator_m

program boxFlowField
    use conditionValue_m
    use unstructuredGrid_mod
    use terminalControler_m
    use boxOperator_m
    implicit none
    integer i_box, nc, num_box
    character(255) caseName
    type(boxResult_t), allocatable :: bResult(:)
    type(UnstructuredGrid) mesh

    call setup_boxArray

    num_box = size(box_array)

    do nc = 1, iargc()
        ! コマンドライン引数を取得
        call getarg(nc, caseName)
        print*, trim(caseName)

        ! mesh = UnstructuredGrid_(condVal%path2FlowFile, condVal%meshFile)
        if(nc==1) then
            mesh = UnstructuredGrid_(trim(caseName)//'/field_0000005125.array', './case1.vtk')
        else
            call mesh%updateWithFlowFieldFile(trim(caseName)//'/field_0000005125.array')
        end if

        allocate(bResult(num_box))

        call set_formatTC('("BoxCellSerch [ #box : ",i6," / ",i6," ]")')

        block
            integer i_cell
            i_cell = 1
            do i_box = 1, num_box
                call print_progress([i_box, num_box])
                ! i_cell = mesh%nearest_cell(box_array(i_box)%center)
                call mesh%search_refCELL(box_array(i_box)%center, i_cell)
                bResult(i_box)%flowVelocity = mesh%CELLs(i_cell)%flowVelocity
            end do

        end block
    
        call output_flowCSV(trim(caseName)//'/BoxFlow.csv', bResult)
        call output_boxVTK(trim(caseName)//'/BoxFlow.vtk', bResult)

        deallocate(bResult)

    end do

end program boxFlowField