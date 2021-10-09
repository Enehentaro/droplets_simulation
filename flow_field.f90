module flow_field
    use unstructuredGrid_mod
    use CUBE_mod
    implicit none

    integer INTERVAL_FLOW                           !気流データ出力間隔

    character PATH_AIR*99, HEAD_AIR*20, FNAME_FMT*30 !気流データへの相対パス,ファイル名接頭文字,ファイル名形式
    character(3) FILE_TYPE  !ファイル形式（VTK,INP）
    integer, private :: FNAME_DIGITS !ファイル名の整数部桁数

    double precision MAX_CDN(3), MIN_CDN(3)         !節点座標の上限および下限

    type reference_cell_t
        integer :: ID = 0, nodeID(3) = 0
    end type reference_cell_t

    contains
    !***********************************************************************
    subroutine set_FILE_TYPE
        integer i

        i = index(FNAME_FMT, '0')
        HEAD_AIR = FNAME_FMT(: i-1)             !ファイル名の接頭文字(最初のゼロの手前まで)
        FNAME_DIGITS = index(FNAME_FMT, '.') - i   !ファイル名の整数部桁数(最初のゼロの位置からドットまでの文字数)

        if(index(FNAME_FMT, '.vtk') > 0) then
            FILE_TYPE = 'VTK'

        else if(index(FNAME_FMT, '.inp') > 0) then
            FILE_TYPE = 'INP'

        else if(index(FNAME_FMT, '.fld') > 0) then
            FILE_TYPE = 'FLD'

        else if(index(FNAME_FMT, '.f') > 0) then
            FILE_TYPE = 'P3D'

        else
            print*, 'FNAME_FMT NG:', FNAME_FMT
            stop
        end if

        print*, 'FILE_TYPE: ', FILE_TYPE, ' ', trim(FNAME_FMT)

    end subroutine set_FILE_TYPE

    character(4) function get_digits_format()

        write(get_digits_format,'("i", i1, ".", i1)') FNAME_DIGITS, FNAME_DIGITS

    end function get_digits_format

    subroutine pre_setting_onFlow
        if(FILE_TYPE=='P3D') then
            call read_faceShape(PATH_AIR)
            call set_faceShape
        else
            call read_nextcell(PATH_AIR)
            call read_boundaries(PATH_AIR)
        end if

    end subroutine

    function get_velocity_flow(cell) result(velocity)
        type(reference_cell_t) :: cell
        double precision velocity(3)

        if(FILE_TYPE == 'P3D') then
            velocity(:) = get_velocity_f(cell%nodeID, cell%ID)
        else
            velocity(:) = VELC(:, cell%ID)
        end if
    end function get_velocity_flow

    function search_ref_cell(X, ref_cel_pre) result(reference_cell)
        double precision, intent(in) :: X(3)
        type(reference_cell_t), intent(in) :: ref_cel_pre
        type(reference_cell_t) reference_cell

        if(FILE_TYPE=='P3D') then
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

        else
            if(ref_cel_pre%ID == 0) then   !参照セルが見つかっていない（＝初期ステップ）
                reference_cell%ID = nearest_cell(X)    
                print*, 'FirstNCN:', reference_cell%ID
        
            else
                reference_cell%ID = nearer_cell(X, ref_cel_pre%ID)
                if (reference_cell%ID == 0) then
                    print*, 'NCN_ERROR:', X(:), reference_cell
                    stop
                end if
    
                if (.not.nearcell_check(X(:), reference_cell%ID)) reference_cell%ID = nearest_cell(X)
        
            end if

        end if

    end function search_ref_cell
    
end module flow_field