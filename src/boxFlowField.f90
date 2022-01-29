program boxFlowField
    use conditionValue_m
    use boxCounter_m
    use caseName_m
    use unstructuredGrid_mod
    implicit none
    integer i_box, num_box, nc_max
    integer, pointer :: nc => nowCase
    character(255) caseName
    type(conditionValue_t) condVal
    ! type(BasicParameter) baseParam
    type(boxCounter), allocatable :: box_array(:)

    type boxResult_t
        real flowVelocity(3)
    end type
    type(boxResult_t), allocatable :: bResult(:)

    type(UnstructuredGrid) mesh

    call case_check(num_case = nc_max)
    !print*, 'caseName = ?'
    !read(5, *) caseName

    do nc = 1, nc_max
        caseName = get_caseName()
        call condVal%read(trim(caseName))
        ! baseParam = BasicParameter_(condVal%dt, condVal%L, condVal%U)
    
        box_array = get_box_array(trim(caseName), 0)
    
        num_box = size(box_array)

        ! mesh = UnstructuredGrid_(condVal%path2FlowFile, condVal%meshFile)
        mesh = UnstructuredGrid_('/home/master/droplet/office1_960_246/honkeisan_haiki_a/field_0000005125.array', 'case1.vtk')

        allocate(bResult(num_box))

        block
            integer i_cell

            do i_box = 1, num_box
                print*, i_box,'/',num_box
                i_cell = mesh%nearest_cell(box_array(i_box)%center)
                bResult(i_box)%flowVelocity = mesh%CELLs(i_cell)%flowVelocity
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
        use vtkMesh_operator_m
        type(vtkMesh) boxMesh
        integer i, j, k
        real, parameter :: trans(3,8) = reshape([ &
                                            0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, &
                                            0.0,0.0,1.0, 1.0,0.0,1.0, 0.0,1.0,1.0, 1.0,1.0,1.0], shape(trans))
        real velArray(3, size(box_array))
                                        
        call boxMesh%allocation_node(num_box*8)
        call boxMesh%allocation_cell(num_box)

        do i = 1, num_box

            do j = 1, 8
                k = 8*(i-1) + j - 1
                boxMesh%node_array(k)%coordinate(:) = box_array(i)%min_cdn(:) + box_array(i)%width(:)*trans(:,j)
            end do

            boxMesh%cell_array(i-1)%nodeID = [(8*(i-1) + j - 1, j = 1, 8)]
            boxMesh%cell_array(i-1)%n_TYPE = 11

        end do

        do i = 1, num_box

            velArray(:,i) = bResult(i)%flowVelocity

        end do

        call boxMesh%output(trim(caseName)//'/BoxFlow.vtk', cellVector=velArray, vectorName='VEL')

    end subroutine

end program boxFlowField