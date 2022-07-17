program boxFlowField
    use conditionValue_m
    use boxCounter_m
    use unstructuredGrid_m
    use terminalControler_m
    implicit none
    integer i_box, num_box, nc
    character(255) caseName
    type(boxCounter), allocatable :: box_array(:)

    type boxResult_t
        real flowVelocity(3)
    end type
    type(boxResult_t), allocatable :: bResult(:)

    type(FlowFieldUnstructuredGrid) mesh

    box_array = get_box_array('.', 0)
    num_box = size(box_array)

    do nc = 1, iargc()
        ! コマンドライン引数を取得
        call getarg(nc, caseName)
        print*, trim(caseName)

        ! mesh = FlowFieldUnstructuredGrid_(condVal%path2FlowFile, condVal%meshFile)
        if(nc==1) then
            mesh = FlowFieldUnstructuredGrid_(trim(caseName)//'/field_0000005125.array', './case1.vtk')
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
                bResult(i_box)%flowVelocity = mesh%get_flowVelocityInCELL(i_cell)
            end do

        end block
    
        call output_countCSV
        call output_boxVTK

        deallocate(bResult)

    end do

    contains

    subroutine output_countCSV
        integer n_unit, i
        character(:), allocatable :: csvFName

        csvFName = trim(caseName)//'/BoxFlow.csv'
        print*, 'output: ', csvFName

        open(newunit=n_unit, file=csvFName, status='replace')
        
            write(n_unit, '("x,y,z,u,v,w")')
            
            do i = 1, size(box_array)
                write(n_unit,'(*(g0:,","))') box_array(i)%center, bResult(i)%flowVelocity
            end do

        close(n_unit)

    end subroutine

    subroutine output_boxVTK
        use VTK_operator_m
        type(UnstructuredGrid_inVTK) boxMesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
        real velArray(3, size(box_array))
        real, allocatable :: xyz(:,:)
        integer, allocatable :: vertices(:,:), types(:)
                                        
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

        boxMesh = UnstructuredGrid_inVTK_(xyz, vertices, types)

        call boxMesh%output(trim(caseName)//'/BoxFlow.vtk', cellVector=velArray, vectorName='VEL')

    end subroutine

end program boxFlowField