!CUBE格子上の流速場をVTK非構造格子に載せるプログラム。
!非構造格子上の各格子に対して、CUBE格子上の最近傍節点を探し、対応付けを行う。
!対応する各節点における流速を配列にして、そのままバイナリファイル出力を行う。
program CUBE2USG
    use plot3d_operator
    use VTK_operator_m
    use simpleFile_reader
    use array_m
    implicit none
    character(100) F_fname, USG_fname, casefname
    character(50),parameter :: filename = 'name.txt'
    character(50), allocatable :: field_name(:), caseName(:)
    character(20), parameter :: CorrespondenceFName = 'vtkCell2cubeNode.txt'
    integer fileID, num_record, num_cell, nc, nc_max
    real, allocatable :: velocity(:,:)
    type(UnstructuredGrid_inVTK) USG
    type(Plot3dMesh) cubeMesh
    type(plot3dNodeInfo), allocatable :: vtkCell2cubeNode(:)

    call read_textRecord(filename, field_name)
    num_record = size(field_name)

    print *, 'Case Name ?'
    read(5,'(A)') casefname
    call read_textRecord(casefname, caseName)
    nc_max = size(caseName)
    
    print *, 'UnstructuredGRID_FileName ?'
    read(5,*) USG_fname

    call USG%read(USG_fname)

    num_cell = USG%get_numCell()

    cubeMesh = read_plot3d_multigrid('mesh.g')   !Gファイル読み込み
    
    do nc = 1, nc_max

        do fileID = 1, num_record
            F_fname = trim(caseName(nc))//'/output/'//trim(field_name(fileID))

            call cubeMesh%read_plot3d_function(F_fname)   !Fファイル読み込み

            if (.not.allocated(vtkCell2cubeNode)) call solve_correspondence

            if (.not.allocated(velocity)) allocate(velocity(3, num_cell))

            block
                integer cellID
                character(:), allocatable :: fname_base

                do cellID = 1, num_cell
                    velocity(:, cellID) = cubeMesh%get_velocity(vtkCell2cubeNode(cellID))
                end do

                fname_base = F_fname(:len_trim(F_fname)-2)

                !流速場配列をバイナリ出力
                call output_2dArray_asBinary(fname=fname_base//'.array', array=velocity)

                !確認用に、ひとつだけVTKファイル出力
                if(fileID==1) call USG%output(fname_base//'.vtk', cellVector=velocity, vectorName='Velocity')

            end block
        
        end do

    end do

    contains

    !対応する節点情報をアスキーファイルで出力するサブルーチン
    subroutine output_nodeInfo
        integer n_unit, i

        print*, 'output: ', CorrespondenceFName

        open(newunit=n_unit, file=CorrespondenceFName, status='replace', action='write')

            write(n_unit, '("#cube ", i0)') cubeMesh%get_numCube()
            write(n_unit, '("cubeshape", 3(x, i0))') cubeMesh%get_cubeShape()

            write(n_unit, '("#usgcell ", i0)') num_cell
            do i = 1, num_cell
                write(n_unit, '(i0, x, 3(x, i0))') vtkCell2cubeNode(i)%cubeID, vtkCell2cubeNode(i)%nodeID(:)
            end do

        close(n_unit)

    end subroutine

    !非構造格子に対応する節点情報を探すサブルーチン
    subroutine search_nodeInfo
        use timeKeeper_m
        use terminalControler_m
        use unstructuredElement_m
        integer i
        real progress_percent, estimation, speed
        type(TimeKeeper) tk
        real, allocatable :: cellCenter(:,:)

        tk = TimeKeeper_()

        print*, "START : CUBENODE SEARCH"
        call set_formatTC('("SEARCH vtkcell2cubenode [ ",f6.2," % ] ", f8.1, " sec is left.")')

        allocate(vtkCell2cubeNode(num_cell))
        cellCenter = get_cellCenters(USG%node_array, USG%cell_array)
        do i = 1, num_cell

            progress_percent = real(i*100) / real(num_cell)
            speed = real(i) / tk%erapsedTime()
            estimation = real(num_cell - i) / (speed + 1.e-9)
            call print_progress([progress_percent, estimation])

            vtkCell2cubeNode(i) = cubeMesh%nearestNodeInfo(cellCenter(:,i))
            
        end do

    end subroutine

    !節点情報対応付けファイルを読み込むサブルーチン
    subroutine read_nodeInfo(success)
        use array_m
        logical, intent(out) :: success
        integer n_unit, i, num_cell_, num_cube, cubeShape(3)
        character dummy*10

        print*, 'read: ', CorrespondenceFName

        open(newunit=n_unit, file=CorrespondenceFName, status='old', action='read')

            read(n_unit, *) dummy, num_cube
            read(n_unit, *) dummy, cubeShape

            read(n_unit, *) dummy, num_cell_

            if((num_cell_/=num_cell).or.(num_cube/=cubeMesh%get_numCube()).or.&
                .not.all(cubeShape==cubeMesh%get_cubeShape())) then

                print*, 'SizeERROR:', num_cell_, num_cell, num_cube, cubeMesh%get_numCube()
                success = .false.
                return
                
            end if

            allocate(vtkCell2cubeNode(num_cell_))
            do i = 1, num_cell_
                read(n_unit, *) vtkCell2cubeNode(i)%cubeID, vtkCell2cubeNode(i)%nodeID(:)
            end do

        close(n_unit)

        success = .true.

    end subroutine

    !格子と節点の対応付けを解決するサブルーチン
    subroutine solve_correspondence
        logical existance, success

        inquire(file=CorrespondenceFName, exist=existance)

        if(existance) then  !対応付けファイルが存在すれば読み込む
            call read_nodeInfo(success)
        else
            success = .false.   !なければ失敗
        end if

        if(.not.success) then   !対応付けに失敗すれば改めて最近傍節点の探索を行う
            call search_nodeInfo
            call output_nodeInfo
        end if

    end subroutine
    
end program CUBE2USG