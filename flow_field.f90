module flow_field
    use unstructuredGrid_mod
    use CUBE_mod
    implicit none

    integer INTERVAL_FLOW                           !気流データ出力間隔

    character PATH_AIR*99, HEAD_AIR*20, FNAME_FMT*30 !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    integer, private :: FNAME_DIGITS !ファイル名の整数部桁数
    logical unstructuredGrid

    double precision MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限

    type reference_cell_t
        integer :: ID = 0, nodeID(3) = 0
    end type reference_cell_t

    contains
    !***********************************************************************
    subroutine check_FILE_GRID
        integer i
        character(3) type_name

        unstructuredGrid = .true.

        i = index(FNAME_FMT, '0')
        HEAD_AIR = FNAME_FMT(: i-1)             !ファイル名の接頭文字(最初のゼロの手前まで)
        FNAME_DIGITS = index(FNAME_FMT, '.') - i   !ファイル名の整数部桁数(最初のゼロの位置からドットまでの文字数)
        
        if(index(FNAME_FMT, '.f') > 0) then
            type_name = 'P3D'
            unstructuredGrid = .false.
            print*, 'FILE_GRID : CUBE', trim(FNAME_FMT)

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
        use adjacent_information
        logical success

        if(unstructuredGrid) then
            call read_adjacency(PATH_AIR, success)
            if(success) then
                call read_boundaries(PATH_AIR)

            else
                call solve_adjacentInformation
                call output_boundaries(PATH_AIR)
                call output_adjacency(PATH_AIR)

            end if
            call boundary_setting

        else
            call read_faceShape(PATH_AIR)
            call set_faceShape
        end if

    end subroutine

    subroutine read_flow_data(FNUM)
        integer, intent(in) :: FNUM
        character(99) FNAME
        character(:),allocatable :: digits_fmt

        digits_fmt = get_digits_format()

        if(unstructuredGrid) then
            if (INTERVAL_FLOW == -1) then !定常解析
                FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
                call read_unstructuredGrid(FNAME)
            else
                FNAME = trim(PATH_AIR)//trim(HEAD_AIR)
                call read_unstructuredGrid(FNAME, digits_fmt, FNUM)
            end if
                       
            MAX_CDN(1) = maxval(CDN(1,:))
            MAX_CDN(2) = maxval(CDN(2,:))
            MAX_CDN(3) = maxval(CDN(3,:))
            print*, 'MAX_coordinates=', MAX_CDN(:)
                
            MIN_CDN(1) = minval(CDN(1,:))
            MIN_CDN(2) = minval(CDN(2,:))
            MIN_CDN(3) = minval(CDN(3,:))
            print*, 'MIN_coordinates=', MIN_CDN(:)
                
            call set_gravity_center
                   
            call boundary_setting

        else

            if (INTERVAL_FLOW == -1) then !定常解析
                FNAME = trim(PATH_AIR)//trim(FNAME_FMT)
            else
                write(FNAME,'("'//trim(PATH_AIR)//trim(HEAD_AIR)//'",'//digits_fmt//',".f")') FNUM

            end if
            call read_CUBE_data(FNAME, trim(PATH_AIR))

            block
                real min_max(6)

                min_max = get_minMax_CUBE()

                MIN_CDN(:) = min_max(1:3)
                MAX_CDN(:) = min_max(4:6)

            end block

        end if
            
    end subroutine read_flow_data

    function search_ref_cell(X, ref_cel_pre) result(reference_cell)
        double precision, intent(in) :: X(3)
        integer, intent(in) :: ref_cel_pre
        integer reference_cell

        if(ref_cel_pre == 0) then   !参照セルが見つかっていない（＝初期ステップ）
            reference_cell = nearest_cell(X)    
            print*, 'FirstNCN:', reference_cell
    
        else
            reference_cell = nearer_cell(X, ref_cel_pre)
            if (reference_cell == 0) then
                print*, 'NCN_ERROR:', X(:), reference_cell
                stop
            end if

            if (.not.nearcell_check(X(:), reference_cell)) reference_cell = nearest_cell(X)
    
        end if

    end function search_ref_cell

    function search_ref_cell_onCUBE(X, ref_cel_pre) result(reference_cell)
        double precision, intent(in) :: X(3)
        type(reference_cell_t), intent(in) :: ref_cel_pre
        type(reference_cell_t) reference_cell

        if(ref_cel_pre%ID == 0) then   !参照セルが見つかっていない（＝初期ステップ）
            reference_cell%ID = get_cube_contains(real(X))    
            reference_cell%nodeID(:) = nearest_node(real(X), reference_cell%ID)
            print*, 'FirstNCN:', reference_cell
    
        else
            reference_cell%ID = ref_cel_pre%ID
            reference_cell%nodeID(:) = nearer_node(real(X), ref_cel_pre%nodeID, ref_cel_pre%ID)

            if (.not.nearNode_check(real(X), reference_cell%nodeID, reference_cell%ID)) then
                reference_cell%ID = get_cube_contains(real(X))    
                reference_cell%nodeID(:) = nearest_node(real(X), reference_cell%ID)
            end if
    
        end if

    end function search_ref_cell_onCUBE
    
end module flow_field