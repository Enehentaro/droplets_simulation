module flow_field
    use unstructuredGrid_mod
    use CUBE_mod
    implicit none

    integer INTERVAL_FLOW                           !気流データ出力間隔
    integer LoopS, LoopF, OFFSET

    character PATH_FlowDIR*99, HEAD_AIR*20, FNAME_FMT*30 !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    integer, private :: FNAME_DIGITS !ファイル名の整数部桁数
    
    logical unstructuredGrid

    real MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限

    type reference_cell_t
        integer :: ID = 0, nodeID(3) = 0
    end type reference_cell_t

    contains
    !***********************************************************************
    subroutine check_FILE_GRID
        integer i, i_

        unstructuredGrid = .true.

        i = index(FNAME_FMT, '0')
        i_ = index(FNAME_FMT, '_', back=.true.)
        if(i_ > i) i = i_ + 1
        HEAD_AIR = FNAME_FMT(: i-1)             !ファイル名の接頭文字(最初のゼロの手前まで)
        FNAME_DIGITS = index(FNAME_FMT, '.') - i   !ファイル名の整数部桁数(最初のゼロの位置からドットまでの文字数)

        i = len_trim(FNAME_FMT)
        if(FNAME_FMT(i-1 : i) == '.f') then
            unstructuredGrid = .false.
            print*, 'FILE_GRID : [CUBE] ', trim(FNAME_FMT)

        else
            call check_FILE_TYPE(FNAME_FMT)

        end if

    end subroutine check_FILE_GRID

    function get_digits_format() result(format)
        character(2*(FNAME_DIGITS/10 +1) + 2) format

        if(FNAME_DIGITS <= 9) then
            write(format,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS
        else
            write(format,'("i", i2, ".", i2)') FNAME_DIGITS, FNAME_DIGITS
        end if

    end function get_digits_format

    subroutine preprocess_onFlowField
        use adhesion_onSTL_m
        use adjacent_information
        logical success

        if(unstructuredGrid) then
            call read_adjacency(PATH_FlowDIR, success)
            if(success) then
                call read_boundaries(PATH_FlowDIR)

            else
                call solve_adjacentInformation
                call output_boundaries(PATH_FlowDIR)
                call output_adjacency(PATH_FlowDIR)

            end if

        else
            call read_faceShape(PATH_FlowDIR)
            call set_faceShape
        end if

    end subroutine preprocess_onFlowField

    subroutine read_steadyFlowData
        character(99) FNAME

        if(unstructuredGrid) then
            FNAME = trim(PATH_FlowDIR)//trim(FNAME_FMT)
            call read_unstructuredGrid(FNAME)
            call set_gravity_center

        else
            FNAME = trim(PATH_FlowDIR)//trim(FNAME_FMT)
            call read_CUBE_data(FNAME, trim(PATH_FlowDIR))

        end if
            
    end subroutine read_steadyFlowData

    subroutine read_unsteadyFlowData(FNUM)
        integer, intent(in) :: FNUM
        character(99) FNAME
        character(:),allocatable :: digits_fmt

        digits_fmt = get_digits_format()

        if(unstructuredGrid) then
            FNAME = trim(PATH_FlowDIR)//trim(HEAD_AIR)
            call read_unstructuredGrid(FNAME, digits_fmt, FNUM)
            call set_gravity_center

        else
            write(FNAME,'("'//trim(PATH_FlowDIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".f")') FNUM
            call read_CUBE_data(FNAME, trim(PATH_FlowDIR))

        end if
            
    end subroutine read_unsteadyFlowData

    subroutine set_MinMaxCDN
        real min_max(6)

        if(unstructuredGrid) then
            MAX_CDN(1) = maxval(NODEs(:)%coordinate(1))
            MAX_CDN(2) = maxval(NODEs(:)%coordinate(2))
            MAX_CDN(3) = maxval(NODEs(:)%coordinate(3))
            print*, 'MAX_coordinates=', MAX_CDN(:)
                
            MIN_CDN(1) = minval(NODEs(:)%coordinate(1))
            MIN_CDN(2) = minval(NODEs(:)%coordinate(2))
            MIN_CDN(3) = minval(NODEs(:)%coordinate(3))
            print*, 'MIN_coordinates=', MIN_CDN(:)

        else
            min_max = get_minMax_CUBE()

            MIN_CDN(:) = min_max(1:3)
            MAX_CDN(:) = min_max(4:6)
        end if
            

    end subroutine set_MinMaxCDN

    subroutine search_refCELL(X, reference_cell)
        real, intent(in) :: X(3)
        integer, intent(inout) :: reference_cell

        if(reference_cell == 0) then   !参照セルが見つかっていない（＝初期ステップ）
            reference_cell = nearest_cell(X)
            print*, 'FirstNCN:', reference_cell
    
        else
            reference_cell = nearer_cell(X, reference_cell)
            if (reference_cell == 0) then
                print*, 'NCN_ERROR:', X(:), reference_cell
                stop
            end if

            if (.not.nearcell_check(X(:), reference_cell)) reference_cell = nearest_cell(X)
    
        end if

    end subroutine search_refCELL

    subroutine search_refCELL_onCUBE(X, reference_cell)
        real, intent(in) :: X(3)
        type(reference_cell_t), intent(inout) :: reference_cell

        if(reference_cell%ID == 0) then   !参照セルが見つかっていない（＝初期ステップ）
            reference_cell%ID = get_cube_contains(X)    
            reference_cell%nodeID(:) = nearest_node(X, reference_cell%ID)
            print*, 'FirstNCN:', reference_cell
    
        else
            reference_cell%nodeID(:) = nearer_node(X, reference_cell%nodeID, reference_cell%ID)
            if (.not.nearNode_check(X, reference_cell%nodeID, reference_cell%ID)) then
                reference_cell%ID = get_cube_contains(X)    
                reference_cell%nodeID(:) = nearest_node(X, reference_cell%ID)
            end if
    
        end if

    end subroutine search_refCELL_onCUBE

    integer function get_FileNumber(step)
        integer, intent(in) :: step

        get_FileNumber = OFFSET
        do while(get_FileNumber + INTERVAL_FLOW < Step)
            get_FileNumber = get_FileNumber + INTERVAL_FLOW
        end do

    end function get_FileNumber

    subroutine deallocation_flow
        if(unstructuredGrid) then
            call deallocation_unstructuredGRID
        else
            call deallocation_CUBE
        end if
    end subroutine deallocation_flow
    
end module flow_field