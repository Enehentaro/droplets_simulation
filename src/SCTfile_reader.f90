module SCT_file_reader_m
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !!package     : SCT_file_reader ver 2.5.2
    !!author      : T.Ikeda, Y.Ida
    !!description :
    !! - SC/TETRA 出力のファイルを読み取り，データを取り出す.  
    !! - 今までのconverterと異なり, 本モジュールで独立して扱えるようになっている. 
    !! - 並列化には対応していない. 
    !!note        :
    !! - セルの節点の並び順はwedge以外はvtkのものと同じ.
    !! - SC/TETRAでは，セル番号および節点番号は0スタートなので, fortran運用のためインデックス+1. 
    !! - 頂点配列(NDNO)はセルタイプ毎に並んでいない. 
    !!Others      :
    !! - cell2verticesのrank 1には最大で8つ(hexahedron), face2verticesには4つの値が入るが, 
    !!   値が入っていない箇所は全て-1に統一されている. 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none
    private

    integer(4) LCORD
        !! 座標番号(=0固定): 多分不要

    integer,parameter :: LengthOfSubRecord = 4194304
        !! サブレコード1行あたりの長さの上限. 必要はないがメモとして

    !面の定義. [頂点数, ローカル節点番号の並び].ローカル節点番号は1~最大8.  
    integer, private, target :: SCTTetraFaces_(3+1,4)   = reshape( &
                                                    [3, 4,3,2, &
                                                     3, 4,1,3, &
                                                     3, 4,2,1, &
                                                     3, 1,2,3], shape(SCTTetraFaces_))
    integer, private, target :: SCTPyramidFaces_(4+1,5) = reshape( &
                                                    [3, 1,5,2,-1, &
                                                     3, 2,5,3,-1, &
                                                     3, 3,5,4,-1, &
                                                     3, 4,5,1,-1, &
                                                     4, 1,2,3,4], shape(SCTPyramidFaces_))
    integer, private, target :: SCTPrismFaces_(4+1,5)   = reshape( &
                                                    [4, 1,4,5,2, &
                                                     4, 2,5,6,3, &
                                                     4, 3,6,4,1, &
                                                     3, 1,2,3,-1, &
                                                     3, 6,5,4,-1], shape(SCTPrismFaces_))        
    integer, private, target :: SCTHexahedronFaces_(4+1,6)   = reshape( &
                                                    [4, 1,5,6,2, &
                                                     4, 2,6,7,3, &
                                                     4, 3,7,8,4, &
                                                     4, 4,8,5,1, &
                                                     4, 1,2,3,4, &
                                                     4, 8,7,6,5], shape(SCTHexahedronFaces_))

    real(8),parameter,public :: MissingValueSize = 1.0d20
        !! SC/TETRAで規定された欠測値の大きさ. 

    character(10),parameter :: LS_Scalar_HeadName = "LS_Scalar:"
    character(10),parameter :: LS_Vector_HeadName = "LS_Vector:"
    character(32),parameter :: OverLapEndLabel = "OverlapEnd"

    integer,parameter :: LS_Scalar_DataSize = 30
        !!.fldに含まれるスカラーデータの個数上限. 
    integer,parameter :: LS_Vector_DataSize = 30
        !!.fldに含まれるベクトルデータの個数上限. 

    type,private :: LS_Scalar_t
        character(:),allocatable :: name
        character(:),allocatable :: abbreviated_name
        integer(4) ndata
        real(8),allocatable :: data(:)
    end type

    type,private :: LS_Vector_t
        character(:),allocatable :: name
        character(:),allocatable :: abbreviated_name
        integer(4) ndata
        real(8),allocatable :: x(:), y(:), z(:)
    end type
    type,private :: sctregion_t
        integer(4) NE
        character(:),allocatable :: LRGN
        integer(4),allocatable :: IE(:)
            !!領域を構成する要素番号
        integer(4),allocatable :: IFA(:)
            !!領域を構成するローカル面番号. 
        contains
        procedure,private :: get_data => get_sctregion_data_
        procedure,private :: extract_region_surface_
    end type sctregion_t

    type sct_data_name_list_t
        !! スカラーorベクトルデータの名前だけを取り出す. 構造体の配列にして使用する.
        !! region用に使うことも出来る. その場合abbreviatedは使わない.
        !! 異なる文字長の配列が実装できなかったのでこれで代用する.  
        character(:),allocatable :: name
        character(:),allocatable :: abbreviated_name
    end type

    type sct_grid_t
        !! SC/TETRA メッシュクラス. 
        !! 必要最低限の変数のみ保持. 変数名はフォーマットに準拠.  
        !! メッシュそのものを取り扱うのでメモリ圧迫する可能性大. 
        !! ソルバ内で使う場合はサブルーチンのローカル変数として扱う方が無難(自動開放されるはず)
        logical,private :: is_FLD
        logical,private :: includes_topo
        ! integer,private :: error_code

        integer(4),private ::  NELEM = 0
            !!要素数
        integer(4),private ::  NTTE = 0
            !!1要素あたりの節点数の合計

        integer(4),allocatable,private :: IETYP(:)
            !! 要素タイプ. 
            !! 34:tetrahedron, 35:pyramid, 36:wedge, 38:hexahedron
        integer(4),allocatable,private :: NDNO(:)
            !! 節点番号

        integer(4),allocatable,private :: GRP(:)
            !! グループ番号

        integer(4),allocatable,private :: MAT(:)
            !! 物性番号
            !! 多分使えないけど念のため...

        integer(4),private :: NNODS = 0
            !! 節点の総数
        real(8),allocatable,private :: CDN_X(:), CDN_Y(:), CDN_Z(:)
            !! 節点座標

        integer(4),private ::  NDATA = 0
            !! 場の変数のデータ数(=NNODS)

        ! SctRegion用変数
        integer(4),private :: NRGN = 0
            !! NRGN: 領域の個数, 配列regionの次元
        integer(4),private :: NLEN = 0
            !! NLEN: 領域名LRGNのバイト数 

        type(sctregion_t),allocatable,private :: region(:)
        
        !test 2021/10/12
        !連結リストの方が良い?
        integer,private :: scalar_data_count = 0
        integer,private :: vector_data_count = 0
        type(LS_Scalar_t),allocatable,private :: scalars(:)
        type(LS_Vector_t),allocatable,private :: vectors(:)

        integer,private :: tetra_count_ = 0
        integer,private :: pyramid_count_ = 0
        integer,private :: wedge_count_ = 0
        integer,private :: hexa_count_ = 0

        contains

        procedure,public :: includes_topology
        procedure,public :: is_fld_file

        procedure,public :: print_self

        procedure,public :: read_SCT_file

        procedure,public :: extract_original_cell_vertices
        procedure,public :: extract_cell_vertices
        procedure,public :: extract_ordered_cell_vertices

        procedure,public :: get_2d_array_of_point_coords
        ! procedure get_2d_array_of_point_velocity !データをsearch_vector_dataで探す方針なので削除. 

        procedure,public :: get_cell_types

        procedure,public :: get_element_count
        procedure,public :: get_vertex_count

        procedure,public :: get_tetrahedron_count
        procedure,public :: get_wedge_count
        procedure,public :: get_pyramid_count
        procedure,public :: get_hexahedron_count

        procedure,public :: get_region_count
        procedure,public :: get_region_namelist
        procedure,public :: extract_face2vertices_on_region

        procedure,public :: search_scalar_data
        procedure,public :: search_vector_data
        procedure get_data_titles

        final destructor

    end type

    interface get_data_array_
        module procedure get_data_array_int32_
        module procedure get_data_array_float64_
    end interface get_data_array_


    public sct_grid_t, sct_data_name_list_t, get_cell_data_from_cellvertices
contains

    !==========================================================================
    ! PUBLIC : subroutine
    !==========================================================================


    subroutine get_cell_data_from_cellvertices(cell_data, cell2vertices, point_data)
        !! 各セル毎の頂点配列に関連する節点中心データからセル中心データを構築する. 
        !! 値はセルを構成する節点データの算術平均として計算する. 
        
        implicit none
        integer,intent(inout) :: cell_data(:)
            !!出力されるセル中心データ. 
        integer,intent(in) :: cell2vertices(:,:)
            !!頂点配列. 
        integer,intent(in) :: point_data(:)
            !!任意の節点データ. 

        integer cell, node_count
        
        node_count = ubound(cell2vertices,dim=1)

        do cell = 1, size(cell2vertices,dim=2)
            cell_data(cell) = sum(point_data(cell2vertices(1:node_count,cell))) / node_count
        end do
        
    end subroutine

    !==========================================================================
    ! type bounded procedure : sct_grid_t
    !==========================================================================
    logical function includes_topology(this)
        !! 格子ファイルがトポロジを含むか. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        includes_topology = this%includes_topo
    
    end function

    logical function is_fld_file(this)
        !!ファイルがFLDか.
        implicit none
        class(sct_grid_t),intent(in) :: this

        is_fld_file = this%is_FLD

    end function

    subroutine print_self(this,unit)
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,intent(in) :: unit
        integer i, nrg

        write(unit,"('sct grid information ')",advance="no") 
        if(this%includes_topo) then
            write(unit, "(A)") "from "//merge("fld file", "pre file",this%is_FLD)
            write(unit, "(' cell      :: ', i0)") this%NELEM
            write(unit, "(' - tetra   :: ', i0)") this%tetra_count_
            write(unit, "(' - pyramid :: ', i0)") this%pyramid_count_
            write(unit, "(' - wedge   :: ', i0)") this%wedge_count_
            write(unit, "(' - hexa    :: ', i0)") this%hexa_count_
            write(unit, "(' node      :: ', i0)") this%NNODS

            if(allocated(this%region)) then
                write(unit, "(' region info ')")
                do nrg = 1, this%NRGN
                    write(unit, "(2x,A,1x,A)") trim(this%region(nrg)%LRGN), merge("(vol )", "(surf)", all(this%region(nrg)%IFA==0))
                end do
            endif

        else
            write(unit, "('topology information not found.')")
        endif

        if(this%is_FLD) then
            write(unit, "(' scalar data  :: ', i0)") this%scalar_data_count
            write(unit, "(A,' = ',A)") (this%scalars(i)%abbreviated_name, this%scalars(i)%name, i = 1, this%scalar_data_count)
            write(unit, "(' vector data  :: ', i0)") this%vector_data_count
            write(unit, "(A,' = ',A)") (this%vectors(i)%abbreviated_name, this%vectors(i)%name, i = 1, this%vector_data_count)
        endif

    end subroutine

    subroutine read_SCT_file(this, filename)
        !! SCTファイルを開き，データを取得する. 事実上のコンストラクタ. 
        !! 既に別のファイルを開いていた場合，そのデータを破棄して開く. 
        implicit none
        class(sct_grid_t),intent(inout) :: this
        character(*),intent(in) :: filename
        
        ! character(256) :: errmsg

        integer unit

        !割り付けされている物があれば解放する. 
        call destructor(this)


        if(.not. open_binary_sequential_(unit, filename)) then
            return 
        endif

        select case(get_extension_(filename))

        case(".pre")
            ! メッシュファイルの場合
            this%is_FLD = .false.
            this%includes_topo = .true.
            call readPRE_Header_data_(unit)
            call readPRE_Main_data_(unit, this%NELEM, this%IETYP, this%NTTE, this%NDNO, this%GRP, this%MAT, this%NNODS, &
                                   this%CDN_X, this%CDN_Y, this%CDN_Z, this%region, this%NRGN, this%NLEN)

        case(".fld")
            ! fldファイルの場合
            ! fldの初期ファイルでないとトポロジー情報が入手できないらしい. 
            this%is_FLD = .true.
            call readFLD_Header_data_(unit)
            !call readFLD_Main_data_(unit, this%NELEM, this%IETYP, this%NTTE, this%NDNO, this%MAT, this%NNODS, &
            !                       this%CDN_X, this%CDN_Y, this%CDN_Z, this%NDATA, this%VEL_X, this%VEL_Y, this%VEL_Z, this%PRES)
        
            call readFLD_Main_data_2(unit, this%NELEM, this%IETYP, this%NTTE, this%NDNO, this%MAT, this%NNODS, &
                                     this%CDN_X, this%CDN_Y, this%CDN_Z, this%scalars, this%vectors,           &
                                     this%scalar_data_count, this%vector_data_count, this%includes_topo)
        case default
            print "(A,' is not supported. STOP')", trim(filename)
            return
        end  select

        !各セル数. 
        if ( this%includes_topo  ) then
            this%tetra_count_ = count(this%IETYP==34)
            this%pyramid_count_ = count(this%IETYP==35)
            this%wedge_count_ = count(this%IETYP==36)
            this%hexa_count_ = count(this%IETYP==38)
        end if
    
        close(unit)

        ! if(this%error_code == -3) stop "invalid data format. "
    end subroutine 

    subroutine extract_original_cell_vertices(this, cell2vertices)
        !!Sc/Tetraで出力されたセル-頂点関係の配列をそのまま出力する. 
        !!セルの種類毎に並んでいないのが特徴.  頂点番号は1から始まる.
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,allocatable,intent(inout) :: cell2vertices(:,:)
            !!頂点配列. 1st arg: vertex count, 2nd arg: cell number

        integer II
        integer offset
        integer nc
        integer KK_beg, KK_end

        if(.not. this%includes_topo) return

        if(.not. allocated(cell2vertices)) allocate(cell2vertices(1:8,1:this%NELEM),source = -1)

        nc = 1
        KK_beg = 1
        do II = 1, this%NELEM

            offset = this%IETYP(II) - 30
            KK_end = KK_beg + offset - 1

            cell2vertices(1:offset,nc) =this%NDNO(KK_beg:KK_end)
            nc = nc + 1
           
            KK_beg = KK_end + 1

        end do


    end subroutine


    subroutine extract_cell_vertices(this, tetras, pyramids, wedges, hexas)
        !!afdet solver との互換性のため, セルタイプごとの頂点配列を出力する. 
        !!頂点配列にはセル毎の頂点のインデックスが格納される. 頂点番号は1から始まる. 
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,allocatable,intent(inout),optional :: tetras(:,:), pyramids(:,:), wedges(:,:), hexas(:,:)

        if(.not. this%includes_topo) return

        if(present(tetras)) call extract_primitives_(this, tetras, "tetra")
        if(present(pyramids)) call extract_primitives_(this, pyramids, "pyramid")
        if(present(wedges)) call extract_primitives_(this, wedges, "wedge")
        if(present(hexas)) call extract_primitives_(this, hexas, "hexa")

    end subroutine


    subroutine extract_ordered_cell_vertices(this, cell2vertices)
        !!セルの種類毎に並んだ格子全体の頂点配列を作成する. 
        !!セルはtetra→pyramid→wedge→hexaの順に並べられる. 頂点番号は1から始まる.
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,allocatable,intent(inout) :: cell2vertices(:,:)

        integer,allocatable :: tmp(:,:)
        integer :: celltypes(4) = [34,35,36,38]
        character(5),dimension(4) :: typename = ["tetra", "pyram", "wedge", "hexah"] 
            !gfortranだと文字数が同じでないとerrorになるため無理矢理揃えた.  
            !ifortなら問題ないのだが...
        integer i, cell_beg, cell_end
       
        if(.not. this%includes_topo) return
        if(.not. allocated(cell2vertices)) allocate(cell2vertices(1:8,1:this%NELEM))

        cell_beg = 1
        do i = 1, 4
            if(count(this%IETYP==celltypes(i)) > 0) then

                call extract_primitives_(this, tmp, typename(i))

                cell_end = cell_beg + size(tmp, dim=2) - 1

                cell2vertices(1:size(tmp, dim=1), cell_beg:cell_end) = tmp(:,:)

                cell_beg = cell_end + 1

                deallocate(tmp)
            
            end if
        end do    
    
    end subroutine


    subroutine get_2d_array_of_point_coords(this, xyz)
        !! 節点座標を2次元配列で出力する. 
        implicit none
        class(sct_grid_t),intent(in) :: this
        real(8),allocatable,intent(inout) :: xyz(:,:)

        if(.not. this%includes_topo) return
        call packing_vector_into_2Darray_(xyz, this%CDN_X, this%CDN_Y, this%CDN_Z)
    
    end subroutine 


    ! subroutine get_2d_array_of_point_velocity(this, velocity)
    !     !! 節点流速データを2次元配列で出力する. 
    !     implicit none
    !     class(sct_grid_t),intent(in) :: this
    !     real(8),allocatable,intent(inout) :: velocity(:,:)

    !     if(this%error_code /= 0) return
    !     call packing_vector_into_2Darray_(velocity, this%VEL_X, this%VEL_Y, this%VEL_Z)
    
    ! end subroutine 

    subroutine get_cell_types(this, celltypes,conversion)
        !!セルタイプ配列を出力する. 
        !!extract_original_cell_verticesで出力したセル-節点配列に対して有効. 
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,allocatable,intent(inout) :: celltypes(:)
            !!セルタイプ配列
        character(*),intent(in),optional :: conversion
            !!セルタイプ番号をvtk, xdmfいずれかに変換する. 

        if(.not. this%includes_topo) return

        allocate(celltypes, source = this%IETYP)

        if(.not.present(conversion)) return

        select case(trim(conversion))
            case("vtk","VTK")
                where(celltypes == 34)
                    celltypes = 10
                elsewhere(celltypes == 35)
                    celltypes = 14
                elsewhere(celltypes == 36)
                    celltypes = 13
                elsewhere(celltypes == 38)
                    celltypes = 12
                endwhere
            case("xdmf","XDMF")
                where(celltypes == 34)
                    celltypes = 6
                elsewhere(celltypes == 35)
                    celltypes = 7
                elsewhere(celltypes == 36)
                    celltypes = 8
                elsewhere(celltypes == 38)
                    celltypes = 9
                endwhere
            case default
                print*, trim(conversion), " is not implemented. "
                return
        end select
    
    end subroutine 



    integer function get_element_count(this)
        !!要素数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_element_count = this%NELEM
    
    end function


    integer function get_vertex_count(this)
        !!節点数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_vertex_count = this%NNODS

    end function

    integer function get_tetrahedron_count(this)
        !!格子に含まれるテトラ格子数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_tetrahedron_count = this%tetra_count_
    end function

    integer function get_pyramid_count(this)
        !!格子に含まれるピラミッド格子数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_pyramid_count = this%pyramid_count_
    end function

    integer function get_wedge_count(this)
        !!格子に含まれるプリズム格子数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_wedge_count = this%wedge_count_
    end function

    integer function get_hexahedron_count(this)
        !!格子に含まれるヘキサ格子数を取得する. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        if(this%includes_topo) get_hexahedron_count = this%hexa_count_
    end function

    integer function get_region_count(this)
        !!領域の個数. 
        implicit none
        class(sct_grid_t),intent(in) :: this

        get_region_count = this%NRGN

    end function

    subroutine get_region_namelist(this, name_list)
        implicit none
        class(sct_grid_t),intent(in) :: this
        type(sct_data_name_list_t),allocatable,intent(inout) :: name_list(:)

        integer nrg

        allocate(name_list(1:this%NRGN))
        do nrg = 1, this%NRGN
            name_list(nrg)%name = this%region(nrg)%LRGN
        enddo

    end subroutine


    subroutine extract_face2vertices_on_region(this, region_num, cell2vertices, face2vertices)
        !!任意のregionを構成する頂点配列を取得する. 体積領域は無視する. 
        !!cell2verticesはoriginalの物でなければならない. 
        implicit none
        class(sct_grid_t),intent(in) :: this
        integer,intent(in) :: region_num
            !!region番号. 
        integer,allocatable,intent(in) :: cell2vertices(:,:)
            !!並べ替えのされていないセル-頂点配列. 
        integer,allocatable,intent(inout) :: face2vertices(:,:)
            !!regionを構成する面-頂点配列. 

        ! integer nrg
        if(.not. this%includes_topo) return
        if(region_num > this%NRGN) then
            print*, "error(sct_grid_t) :: region_num must be less than ",this%NRGN
            return
        endif

        call this%region(region_num)%extract_region_surface_(this%IETYP, cell2vertices, face2vertices)
    end subroutine



    subroutine search_scalar_data(this, key, scalar)
        !! .fldに含まれるスカラー場データを取得する.
        !! keyにタイトル名を入れて検索する. 該当しない場合含まれるデータ一覧を表示.  
        implicit none
        class(sct_grid_t),intent(in) :: this
        character(*),intent(in) :: key
            !!取り出したいデータのSC/TETRAでの名称. 
        real(8),allocatable,intent(inout) :: scalar(:)

        integer i

        if(.not. this%is_FLD) return !FLDファイルなら場のデータは含まれると判断. 

        do i = 1, size(this%scalars)
            if ( trim(this%scalars(i)%abbreviated_name) == trim(key) ) then
                if ( .not. allocated(scalar) ) then
                    allocate(scalar, source = this%scalars(i)%data)
                else
                    scalar(:) = this%scalars(i)%data(:)
                end if
                return
            end if
        end do

        print "(A,' does not exists.')", key
        print "('Existing Scalar data are listed below:')"
        do i = 1, this%scalar_data_count
            print "(A,1x,':',A)", this%scalars(i)%name, this%scalars(i)%abbreviated_name
        end do

         
    end subroutine

    subroutine search_vector_data(this, key, vector)
        !! .fldに含まれるベクトル場データを取得する.
        !! keyにタイトル名を入れて検索する. 該当しない場合含まれるデータ一覧を表示. 
        implicit none
        class(sct_grid_t),intent(in) :: this
        character(*),intent(in) :: key
            !!取り出したいデータのSC/TETRAでの名称.
        real(8),allocatable,intent(inout) :: vector(:,:)

        integer i

        if(.not. this%is_FLD) return

        do i = 1, size(this%vectors)
            if ( trim(this%vectors(i)%abbreviated_name) == trim(key) ) then
                call packing_vector_into_2Darray_(vector, this%vectors(i)%x, this%vectors(i)%y, this%vectors(i)%z)
                return
            end if
        end do

        print "(A,' does not exists.')", key
        print "('Existing Vector data are listed below:')"
        do i = 1, this%vector_data_count
            print "(A,1x,':',A)", this%vectors(i)%name, this%vectors(i)%abbreviated_name
        end do

    end subroutine

    subroutine get_data_titles(this, titles, data_type)
        !! .fldに含まれるデータのタイトルを取得する. 
        implicit none
        class(sct_grid_t), intent(in) :: this
        type(sct_data_name_list_t),allocatable,intent(inout) :: titles(:)
        character(*), intent(in) :: data_type

        integer i

        if(.not. this%is_FLD) return

        select case(data_type)

            case("scalar")
                allocate(titles(this%scalar_data_count))
                do i = 1, this%scalar_data_count
                    titles(i)%abbreviated_name = this%scalars(i)%abbreviated_name
                    titles(i)%name = this%scalars(i)%name
                end do
            case("vector")
                allocate(titles(this%vector_data_count))
                do i = 1, this%vector_data_count
                    titles(i)%abbreviated_name = this%vectors(i)%abbreviated_name
                    titles(i)%name = this%vectors(i)%name
                end do
            case default
        
        end select

    end subroutine

    subroutine destructor(this)
        implicit none
        type(sct_grid_t),intent(inout) :: this

        ! integer i
    
        if (allocated(this%IETYP)) deallocate(this%IETYP)
        if (allocated(this%MAT)) deallocate(this%MAT)
        if (allocated(this%NDNO)) deallocate(this%NDNO)
        if (allocated(this%GRP)) deallocate(this%GRP)
        if (allocated(this%CDN_X)) deallocate(this%CDN_X)
        if (allocated(this%CDN_Y)) deallocate(this%CDN_Y)
        if (allocated(this%CDN_Z)) deallocate(this%CDN_Z)
        ! if (allocated(this%VEL_X)) deallocate(this%VEL_X)
        ! if (allocated(this%VEL_Y)) deallocate(this%VEL_Y)
        ! if (allocated(this%VEL_Z)) deallocate(this%VEL_Z)
        if (allocated(this%region)) deallocate(this%region)
        if(allocated(this%scalars)) deallocate(this%scalars)
        if(allocated(this%vectors)) deallocate(this%vectors)

        this%NELEM = 0
        this%NNODS = 0
        this%NDATA = 0
        this%NTTE = 0
        this%NLEN = 0
        this%NRGN = 0
        this%tetra_count_ = 0
        this%pyramid_count_ = 0
        this%wedge_count_ = 0
        this%hexa_count_ = 0
        this%scalar_data_count = 0
        this%vector_data_count = 0
    end subroutine

    !==========================================================================
    ! PRIVATE : extraction
    !==========================================================================

    subroutine extract_primitives_(sct_grid, primitive, cell_type)
        !!格子からセルタイプ毎の頂点配列を抜き出す. 
        implicit none
        type(sct_grid_t),intent(in) :: sct_grid
        integer,allocatable,intent(inout) :: primitive(:,:)
            !!セルタイプ毎の頂点配列. 1st arg: local node number, 2nd arg: cell number
            !!頂点インデックスが格納される. インデックスは1スタートになっている. 
        character(*),intent(in) :: cell_type
            !!セルの種類. 
            !!"tetra(or tetrahedron)", "prizm(or wedge)", "pyramid", "hexa(or hexahedron)"

        integer(4) :: nc, offset
        integer(4) :: II, node_of_cell
        integer(4) :: KK_beg, KK_end
        integer a_err
        character(512) msg

        nc = 1
        
        select case(trim(cell_type))
            case("tetra","tetrahedron")
                node_of_cell = 4
            case("prizm","wedge")
                node_of_cell = 6
            case("pyramid","pyram")
                node_of_cell = 5
            case("hexa","hexahedron","hexah")
                node_of_cell = 8
            case default
                print "('Unclassifiable cell type was detected. choose cell type below:')"
                print "('tetra or tetrahedron')"
                print "('prizm or wedge      ')"
                print "('pyramid             ')"
                print "('hexa  or hexahedron ')"
                stop
        end select

        if ( count(sct_grid%IETYP==30+node_of_cell) == 0 ) then
            print "('error: cell_type ', A, ' possibly does not exists on your mesh file...')", cell_type
            ! return
        end if

        !仮に対応するセルタイプが存在しなくても割り付けは出来るが，その場合は以下の処理も自動的にスルーされる
        allocate(primitive(1:node_of_cell,count(sct_grid%IETYP==30+node_of_cell)),source = -1, &
                 stat = a_err, errmsg = msg)

        if ( a_err /= 0 ) then
            print "(i0,1x,A)", a_err, msg
            return
        end if

        KK_beg = 1
        do II = 1, size(sct_grid%IETYP)

            offset = sct_grid%IETYP(II) - 30
            KK_end = KK_beg + offset - 1

            if ( sct_grid%IETYP(II) == 30 + node_of_cell ) then
                primitive(1:node_of_cell,nc) = sct_grid%NDNO(KK_beg:KK_end)
                nc = nc + 1
            end if

            KK_beg = KK_end + 1

        end do

    end subroutine extract_primitives_

    ! subroutine get_cell_data_from_point_data_(cell_data, cell_types, node_numbers, point_data)
    !     !! 節点中心データからセル中心データを取り出す
    !     !! セルの値はセルを構成する節点データの平均として扱われる
    !     !!~~~ CAUTION :: THIS IS DEPRICATED METHOD ~~~
    !     !! REASON) 作ってみたが，あとでセルタイプ毎に並び替えるためこの方法では矛盾する. 
    !     !!         セル中心データへの補間は外部で行うとする方が都合が良い. 
    !     !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !     implicit none
    !     integer,allocatable,intent(inout) :: cell_data(:)
    !     integer,intent(in) :: cell_types(:)
    !     integer,intent(in) :: node_numbers(:)
    !     integer,intent(in) :: point_data(:)

    !     integer cell, node_beg, node_end, node_count

    !     allocate(cell_data(size(cell_types)))

    !     node_beg = 1
    !     do cell = 1, size(cell_types)
    !         node_count = cell_types(cell) - 30
    !         node_end = node_beg + node_count - 1

    !         cell_data(cell) = sum(point_data(node_numbers(node_beg:node_end))) / node_count

    !         node_beg = node_end + 1
            
    !     end do


    ! end subroutine 
    
    !==========================================================================
    ! PRIVATE : utility
    !==========================================================================

    logical function open_binary_sequential_(unit, filename) result(is_opened)
        !! バイナリファイルをシーケンシャル形式で開く. 開けない場合.false.
        integer, intent(inout) :: unit
            !! 装置番号. 
        character(*), intent(in) :: filename
             
        character(256) iomessage
        integer iostatus

        open(newunit = unit, file = filename, iostat = iostatus, iomsg = iomessage, &
             form = 'unformatted', access = 'sequential', status = 'old', convert = 'big_endian')

        if ( iostatus == 0 ) then
            is_opened = .true.
        else
            is_opened = .false.
            print "('error:',i0,1x,A)", iostatus, trim(iomessage)
        end if
    end function

    function get_extension_(filename) result(extension_name)
        !! ファイル名の拡張子を取得する. 
        !! e.g.) hoge.f90なら".f90"が返される. "."が2つ以上ある場合は保証しない. 
        character(*),intent(in) :: filename

        character(:),allocatable :: extension_name

        integer extension_start

        extension_start = index(filename, ".", back = .true.)
        extension_name = filename(extension_start:len_trim(filename))

    end function

    subroutine packing_vector_into_2Darray_(array, x, y, z)
        !! 実数のベクトル配列を2次元配列に詰め直す
        implicit none
        real(8), allocatable, intent(out) :: array(:,:)
        real(8), intent(in) :: x(:), y(:), z(:)

        integer size_of_array

        size_of_array = size(x)
        
        if(size_of_array /= size(y) .or. size_of_array /= size(z)) then
            print "('ERROR(SCTfile_reader, packing_): size of x is different from y or z')"
            return
        end if
        if(.not. allocated(array)) allocate(array(3, size_of_array))

        array(1,:) = x(:)
        array(2,:) = y(:)
        array(3,:) = z(:)
    end subroutine
    
    !==========================================================================
    ! PRIVATE: reading header
    !==========================================================================


    subroutine readPRE_Header_data_(unit)
        !!preファイルのヘッダ部分を読み取る
        implicit none
        integer,intent(in) :: unit

        integer(4)  header_num
        character(8)  header_text
        character(32) title_text
    
        read(unit) header_text
        read(unit) header_num
        ! 序文データの読み取り
        ! おそらくvtk化するためには不要なので
        ! 無限ループでHeaderDataEndを検出した時点でexit
        
        do 
            read(unit) title_text
            if(trim(title_text)=='HeaderDataEnd') exit
            read(unit)
            read(unit)
            read(unit)
        end do
    end subroutine 

    subroutine readFLD_Header_data_(unit)
        implicit none
        integer,intent(in) :: unit
        integer(4)  header_num, NNAMS, n
        character(8)  header_text
        character(32) title_text
    
        read(unit) header_text
        read(unit) header_num
        ! 序文データの読み取り
        ! おそらくvtk化するためには不要なので
        ! 無限ループでHeaderDataEndを検出した時点でexit
        
        do 
            read(unit) title_text
            if(trim(title_text)=='HeaderDataEnd') then
                exit
            elseif((trim(title_text)=='Cycle').or.(trim(title_text)=='Unused')) then
                read(unit)
                read(unit)
    
            elseif(trim(title_text)=='Unit:$TEMP') then
                read(unit)
                read(unit)
                read(unit)
                read(unit)
                read(unit) NNAMS
                do n = 1, NNAMS
                    read(unit)
                end do
                read(unit)
                cycle
            end if
    
            read(unit)
            read(unit)
            read(unit)
    
        end do
    end subroutine


    !==========================================================================
    ! PRIVATE: reading main data
    !==========================================================================


    !> 本文データの読み取り．とにかく全て読み取るようにしている
    subroutine readPRE_Main_data_(unit, NELEM, IETYP, NTTE, NDNO, GRP, MAT, NNODS, CDN_X, CDN_Y, CDN_Z, region, NRGN, NLEN)
        implicit none
        integer,intent(in) :: unit
        integer,intent(inout) :: NELEM
        integer,intent(inout) :: NTTE
        integer,allocatable,intent(inout) :: IETYP(:)
        integer,allocatable,intent(inout) :: NDNO(:)
        integer,allocatable,intent(inout) :: GRP(:)
        integer,allocatable,intent(inout) :: MAT(:)
        integer,intent(inout) :: NNODS
        real(8),allocatable,intent(inout) :: CDN_X(:), CDN_Y(:), CDN_Z(:)
        integer,intent(inout) :: NRGN
        integer,intent(inout) :: NLEN
        type(sctregion_t), allocatable, intent(inout) :: region(:)

        integer(4) nrg
        character(32) main_data_title
        character(32) TITLE

        !> OVerlapStart_nの読み取り
        read(unit) main_data_title
        !本文データの読み取り
        ! irecn == 1 なら タイトル内のサブレコードは1つ
        ! iretn == 1 なら サブレコード内のデータは1つ
        ! ibyte == 4:int32, 8:real64, 1:char(32 or 80 or 1)
        ! 本文データにて ibyte = 1となることはなさそう 

        !> LS_CoordinateSystem
        read(unit) TITLE
        read(unit) !4, 1, 1
        read(unit) !LNX
        call get_data_int32_(unit, LCORD)
        read(unit) !0, 0, 0

        !> LS_Elements
        read(unit) TITLE 
        read(unit) !4, 1, 1
        read(unit) !LNX
        call get_data_int32_(unit, NELEM) 
       
        call get_data_array_(unit, IETYP, NELEM) 
        
        call get_data_int32_(unit, NTTE)
        
        call get_data_array_(unit, NDNO, NTTE)
        !NDNOは0スタートのため, 都合+1する. 
        NDNO(1:NTTE) = NDNO(1:NTTE) + 1  
        read(unit) !0, 0, 0

        !> LS_GrpOfElements
        read(unit) TITLE 
        read(unit) !4, 1, 1
        read(unit) !LNX
        call get_data_int32_(unit, NELEM)
        call get_data_array_(unit, GRP, NELEM)
        read(unit)

        !> LS_MatOfElements
        read(unit) TITLE 
        read(unit) !4, 1, 1
        read(unit) !LNX
        call get_data_int32_(unit, NELEM)
        call get_data_array_(unit, MAT, NELEM)
        read(unit)

        !> LS_Nodes
        read(unit) TITLE 
        read(unit) !4, 1, 1
        read(unit) !LNX
        call get_data_int32_(unit, NNODS)
        call get_data_array_(unit, CDN_X, NNODS)
        call get_data_array_(unit, CDN_Y, NNODS)
        call get_data_array_(unit, CDN_Z, NNODS)
        read(unit)

        !> OverlapEndならそのまま終了
        !> そうで無い場合，任意データLS_SctRegionsが存在しうるが
        !> vtk化のためには不要と思われるので現状は素通り
        !> 追記：2021/05/30　境界面の指定番号なのでIFACE.DATへ利用可能かも
        read(unit) TITLE 
        if(trim(TITLE)=='OverlapEnd') then
            print *, 'SctRegion is not included in this file. '
        else
            print *, 'SctRegions is included.'
            read(unit) !4, 1, 1
            read(unit) !LNX
            call get_data_int32_(unit, NRGN)

            allocate(region(NRGN))

            call get_data_int32_(unit, NLEN)
            do nrg = 1, NRGN
                call region(nrg)%get_data(unit, NLEN)
            end do
            
            read(unit)
        end if

    end subroutine

    subroutine readFLD_Main_data_(unit,NELEM, IETYP, NTTE, NDNO, MAT, &
                                      NNODS, CDN_X, CDN_Y, CDN_Z,    &
                                      NDATA, VEL_X, VEL_Y, VEL_Z,    &
                                      PRES)
        implicit none
        integer,intent(in) :: unit
        integer,intent(inout) :: NELEM
        integer,intent(inout) :: NTTE
        integer,allocatable,intent(inout) :: IETYP(:)
        integer,allocatable,intent(inout) :: NDNO(:)
        integer,allocatable,intent(inout) :: MAT(:)
        integer,intent(inout) :: NNODS
        real(8),allocatable,intent(inout) :: CDN_X(:), CDN_Y(:), CDN_Z(:)
        integer,intent(inout) :: NDATA
        real(8),allocatable,intent(inout) :: VEL_X(:), VEL_Y(:), VEL_Z(:)
        real(8),allocatable,intent(inout) :: PRES(:)

        integer(4) N, NTRY
        character(32) main_data_title
        character(32) TITLE
    
        !> OVerlapStart_nの読み取り
        read(unit) main_data_title
        !本文データの読み取り
        ! irecn == 1 なら タイトル内のサブレコードは1つ
        ! iretn == 1 なら サブレコード内のデータは1つ
        ! ibyte == 4:int32, 8:real64, 1:char(32 or 80 or 1)
        ! 本文データにて ibyte = 1となることはなさそう
    
        do
    
        read(unit) TITLE ; print *, TITLE
        
        select case (trim(TITLE))
            case('LS_CoordinateSystem')
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, LCORD)
                read(unit) !0, 0, 0
    
            case('LS_SurfaceGeometryArray')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NGFAX
                call ignore_data_(unit) !LLEN
                call ignore_data_(unit) !LRGNS
                call ignore_data_(unit) !NBNNS
                call ignore_data_(unit) !IPTYP
                call ignore_data_(unit) !IPMAT
                call ignore_data_(unit) !NTTSS
                call ignore_data_(unit) !NDFA
                read(unit) !0, 0, 0
    
            case('LS_Elements')
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NELEM) 
                call get_data_array_(unit, IETYP, NELEM) 
                call get_data_int32_(unit, NTTE)
                call get_data_array_(unit, NDNO, NTTE)  
                !NDNOは0スタートのため, 都合+1する. 
                NDNO(1:NTTE) = NDNO(1:NTTE) + 1
                read(unit) !0, 0, 0
    
            case('LS_MatOfElements')
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NELEM)
                call get_data_array_(unit, MAT, NELEM)
                read(unit)
    
            case('LS_Nodes')
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NNODS)
                call get_data_array_(unit, CDN_X, NNODS)
                call get_data_array_(unit, CDN_Y, NNODS)
                call get_data_array_(unit, CDN_Z, NNODS)
                read(unit)
    
            case('LS_VolumeGeometryArray')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NVORG
                call ignore_data_(unit) !LLEN
                call ignore_data_(unit) !LRGNS
                call ignore_data_(unit) !NELES
                call ignore_data_(unit) !IELE
                read(unit)
    
            case('LS_RegionName&Type')
                call ignore_data_(unit) !LNX
                call get_data_int32_(unit, NTRY)
                call ignore_data_(unit) !NLEN
                do N = 1, NTRY
                    call ignore_data_(unit) !ITRY
                    call ignore_data_(unit) !MRGN
                end do
                read(unit)
    
            case('LS_SFile')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NLEN
                call ignore_data_(unit) !TEXT
                read(unit)
    
            !ここを切り分ける. scalar dataは要求するデータのみ探索し, 後は捨てる
            !一致するもののみよみとり，あとはignoreする. 
            case('LS_Scalar:PRES')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call get_data_int32_(unit, NDATA) !NDATA
                call get_data_array_(unit, PRES, NDATA) !VAR
                read(unit)
    
            case('LS_Scalar:TEMP')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Scalar:TURK')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
            
            case('LS_Scalar:TEPS')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Scalar:EVIS')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Scalar:YPLS')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Scalar:HTRC')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Scalar:USTR')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('LS_Vector:VEL')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LVCT
                call ignore_data_(unit) !LNAM
                call get_data_int32_(unit, NDATA) !NDATA
                call get_data_array_(unit, VEL_X, NDATA)
                call get_data_array_(unit, VEL_Y, NDATA)
                call get_data_array_(unit, VEL_Z, NDATA)
                read(unit)
    
            case('LS_Vector:HVEC')
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !LVCT
                call ignore_data_(unit) !LNAM
                call ignore_data_(unit) !NDATA
                call ignore_data_(unit) !VAR
                call ignore_data_(unit) !VAR
                call ignore_data_(unit) !VAR
                read(unit)
    
            case('OverlapEnd')
                exit
    
            case default
                print*, 'TITLE_ERROR:', TITLE
                stop
    
        end select
    
        end do
    
    end subroutine

    subroutine readFLD_Main_data_2(unit,NELEM, IETYP, NTTE, NDNO, MAT, &
        NNODS, CDN_X, CDN_Y, CDN_Z,    &
        scalars, vectors, scalar_cnt, vector_cnt, topo_included)
        !!FLDのメインデータ部分を読み取る. 
        !
        !fldが解析スタートファイルでない場合,toporogy情報が含まれない可能性が高い. 
        !したがって, 初期ファイル以外のものを最初に読み取るとトポロジ情報が入手できない．
        !
        implicit none
        integer,intent(in) :: unit
        integer,intent(inout) :: NELEM
        integer,intent(inout) :: NTTE
        integer,allocatable,intent(inout) :: IETYP(:)
        integer,allocatable,intent(inout) :: NDNO(:)
        integer,allocatable,intent(inout) :: MAT(:)
        integer,intent(inout) :: NNODS
        real(8),allocatable,intent(inout) :: CDN_X(:), CDN_Y(:), CDN_Z(:)
        type(LS_Scalar_t),allocatable,intent(inout) :: scalars(:)
        type(LS_Vector_t),allocatable,intent(inout) :: vectors(:)
        integer,intent(inout) :: scalar_cnt
        integer,intent(inout) ::  vector_cnt
        logical,intent(inout) ::  topo_included
            

        integer(4) N, NTRY
        character(32) main_data_title
        character(32) TITLE

        !> OVerlapStart_nの読み取り
        read(unit) main_data_title
        !本文データの読み取り
        ! irecn == 1 なら タイトル内のサブレコードは1つ
        ! iretn == 1 なら サブレコード内のデータは1つ
        ! ibyte == 4:int32, 8:real64, 1:char(32 or 80 or 1)
        ! 本文データにて ibyte = 1となることはなさそう

        topo_included = .false.
        !TITLEは以下の順に並んでいると仮定する. 
        do

            read(unit) TITLE
        
            if( trim(TITLE) == 'LS_CoordinateSystem') then
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, LCORD)
                read(unit) !0, 0, 0

            else if ( trim(TITLE) == 'LS_SurfaceGeometryArray') then
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NGFAX
                call ignore_data_(unit) !LLEN
                call ignore_data_(unit) !LRGNS
                call ignore_data_(unit) !NBNNS
                call ignore_data_(unit) !IPTYP
                call ignore_data_(unit) !IPMAT
                call ignore_data_(unit) !NTTSS
                call ignore_data_(unit) !NDFA
                read(unit) !0, 0, 0

            else if ( trim(TITLE) == 'LS_Elements') then
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NELEM) 
                call get_data_array_(unit, IETYP, NELEM) 
                call get_data_int32_(unit, NTTE)
                call get_data_array_(unit, NDNO, NTTE)  
                !NDNOは0スタートのため, 都合+1する. 
                NDNO(1:NTTE) = NDNO(1:NTTE) + 1
                read(unit) !0, 0, 0

            else if ( trim(TITLE) == 'LS_MatOfElements') then
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NELEM)
                call get_data_array_(unit, MAT, NELEM)
                read(unit)

            else if( trim(TITLE) == 'LS_Nodes') then
                read(unit) !4, 1, 1
                read(unit) !LNX
                call get_data_int32_(unit, NNODS)
                call get_data_array_(unit, CDN_X, NNODS)
                call get_data_array_(unit, CDN_Y, NNODS)
                call get_data_array_(unit, CDN_Z, NNODS)
                read(unit)

                !順番的にこれがトポロジー情報の最後とみる. 
                topo_included = .true.

                !ここでループ抜けすると,LS_Scalar以外の物が前に残っている可能性があるのでNG. 
                ! exit

            else if( trim(TITLE) == 'LS_VolumeGeometryArray') then
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NVORG
                call ignore_data_(unit) !LLEN
                call ignore_data_(unit) !LRGNS
                call ignore_data_(unit) !NELES
                call ignore_data_(unit) !IELE
                read(unit)

            else if ( trim(TITLE) == 'LS_RegionName&Type') then
                call ignore_data_(unit) !LNX
                call get_data_int32_(unit, NTRY)
                call ignore_data_(unit) !NLEN
                do N = 1, NTRY
                    call ignore_data_(unit) !ITRY
                    call ignore_data_(unit) !MRGN
                end do
                read(unit)

            else if( trim(TITLE) == 'LS_SFile') then
                call ignore_data_(unit) !LNX
                call ignore_data_(unit) !NLEN
                call ignore_data_(unit) !TEXT
                read(unit)
                !必須データでないらしい. 確実にここでループ抜け出来るとは限らない. 
                ! exit

            else if( TITLE(1:10) == "LS_Scalar:") then
                !LS_Scalar:以下はループ外で処理するので，1行前に戻ってループ抜けする. 
                backspace(unit)
                exit
                !もしLS_elementsなどを介せずここに来た場合, トポロジは含まれないと判断. 
                 
            else
                print*, "(error) readFLD_main_data_2 :: UnKnown title data was detected. ", trim(TITLE)
                return
            end if

        end do

        !LS_scalar, LS_vector
        if(.not.allocated(scalars)) allocate(scalars(LS_Scalar_DataSize))
        if(.not.allocated(vectors)) allocate(vectors(LS_Vector_DataSize))
        datas_reader:block
            logical EOR
            integer s_cnt, v_cnt, i

            EOR = .false.
            s_cnt = 0
            do while(.not.EOR)
                s_cnt = s_cnt + 1
                call LS_Scalar_reader(scalars(s_cnt), unit, EOR)
            end do
            scalar_cnt = s_cnt - 1

            EOR = .false.
            v_cnt = 0
            do while(.not.EOR)
                v_cnt = v_cnt + 1
                call LS_Vector_reader(vectors(v_cnt), unit, EOR)
            end do
            vector_cnt = v_cnt - 1

            do i = s_cnt, size(scalars)
                scalars(i)%name = "NONE"
            enddo

            do i = v_cnt, size(vectors)
                vectors(i)%name = "NONE"
            enddo

            ! do i = 1, size(scalars)
            !     print *, scalars(i)%name
            !     print *, scalars(i)%abbreviated_name
            ! enddo
    
            ! do i = 1, size(vectors)
            !     print *, vectors(i)%name
            !     print *, scalars(i)%abbreviated_name
            ! enddo
    
        end block datas_reader

    end subroutine

    !==========================================================================
    ! LS_Scalar_t
    !==========================================================================

    subroutine LS_Scalar_reader(ls_scalar, unit, is_ended)
        !!LS_Scalar:で始まる行のデータ読み取り. 
        implicit none
        type(LS_Scalar_t),intent(inout) :: ls_scalar
            !!LS_Scalar_t 構造体
        integer(4),intent(in) :: unit
            !!装置番号
        logical,intent(out) :: is_ended
            !!LS_Scalarの終端かどうか. 
        ! character(:),allocatable,intent(inout) :: LNAME
        ! integer(4),intent(inout) :: NDATA
        ! real(8),allocatable,intent(inout) :: VAR(:)
        
        character(32) title
       
        is_ended = .false.
        read(unit) title
        select case(title(1:10))
            case(LS_Scalar_HeadName)
                ls_scalar%abbreviated_name = title(11:14)
                call ignore_data_(unit)
                call get_data_char_(unit, 32, ls_scalar%name)
                call get_data_int32_(unit, ls_scalar%ndata)
                call get_data_array_(unit, ls_scalar%data, ls_scalar%ndata)
                read(unit)
            case(LS_Vector_HeadName)
                backspace(unit)
                is_ended = .true.
                return
            case(OverLapEndLabel)
                is_ended = .true.
                return
            case default
                print "('something is wrong. ')", trim(title(1:10))
                error stop
        end select
       
    end subroutine

    subroutine LS_Vector_reader(ls_vector, unit, is_ended)
        implicit none
        
        type(LS_Vector_t),intent(inout) :: ls_vector
        integer(4),intent(in) :: unit
        logical,intent(inout) :: is_ended
        ! character(:),allocatable,intent(inout) :: LNAME
        ! integer(4),intent(inout) :: NDATA
        ! real(8),allocatable,intent(inout) :: VAR(:)
        
        character(32) title
       
        is_ended = .false.
        read(unit) title
        select case(title(1:10))
            case(LS_Vector_HeadName)
                ls_vector%abbreviated_name = title(11:14)
                call ignore_data_(unit)
                call ignore_data_(unit)
                call get_data_char_(unit, 32, ls_vector%name)
                call get_data_int32_(unit, ls_vector%ndata)
                call get_data_array_(unit, ls_vector%x, ls_vector%ndata)
                call get_data_array_(unit, ls_vector%y, ls_vector%ndata)
                call get_data_array_(unit, ls_vector%z, ls_vector%ndata)
                read(unit)
            case(OverLapEndLabel)
                is_ended = .true.
                return
            case default
                print "('something is wrong. ')", trim(title(1:10))
                error stop
        end select

    end subroutine
    !==========================================================================
    ! sctregion_t
    !==========================================================================

    !> sctregion_tクラスのメソッド
    subroutine get_sctregion_data_(this, unit, n_length)
        implicit none
        class(sctregion_t),intent(inout) :: this
        integer(4),intent(in) :: unit, n_length

        call get_data_char_(unit, n_length, this%LRGN)
        call get_data_int32_(unit, this%NE)
        call get_data_array_(unit, this%IE, this%NE)
        call get_data_array_(unit, this%IFA, this%NE)

        !要素番号, ローカル面番号は全て1スタートに揃える. 従って体積領域は0番. 
        this%IE(1:this%NE) = this%IE(1:this%NE) + 1
        this%IFA(1:this%NE) = this%IFA(1:this%NE) + 1

    end subroutine 

    subroutine extract_region_surface_(this,celltypes, cell2vertices, face2vert_on_region)
        !! regionを構成する頂点番号の抜き出し. 体積領域はスルー. 
        !! したがってregionの面-頂点関係が取得できる. 
        implicit none
        class(sctregion_t),intent(in) :: this
        integer,intent(in) :: celltypes(:)
        integer,intent(in) :: cell2vertices(:,:)
            !!セル-頂点関係の配列. オリジナルのものでなければならない. 
        integer,allocatable,intent(inout) :: face2vert_on_region(:,:)

        integer i, nv, cell, face
        integer,dimension(:,:),pointer :: SCTElementFacesDef

        if(all(this%IFA == 0)) then
            print*, "this region :"//trim(this%LRGN)//" is not a surface region."
            return
        endif

        if(allocated(face2vert_on_region)) deallocate(face2vert_on_region)
        allocate(face2vert_on_region(1:8,this%NE),source = -1)
        
        do i = 1, this%NE
            cell = this%IE(i)
            face = this%IFA(i)

            select case(celltypes(cell))

                case(34)
                    SCTElementFacesDef => SCTTetraFaces_
                case(35)
                    SCTElementFacesDef => SCTPyramidFaces_
                case(36)
                    SCTElementFacesDef => SCTPrismFaces_
                case(38)
                    SCTElementFacesDef => SCTHexahedronFaces_
                case default
                    print*, "sctregion_t :: not implemented celltype.",celltypes(cell)
                    error stop
            end select

            nv = SCTElementFacesDef(1,face)
            face2vert_on_region(1:nv,i) = cell2vertices(SCTElementFacesDef(2:nv+1,face),cell)
        end do

    end subroutine
    !==========================================================================
    ! data_reader
    !==========================================================================

    !> 整数型データを読みこみ格納
    subroutine get_data_int32_(unit,retval)
        implicit none
        integer(4),intent(in) :: unit
        integer(4),intent(out):: retval
        integer(4) ibyte, iretn, irecn 

        read(unit) ibyte, iretn, irecn
        read(unit) retval 
        
    end subroutine

    !> 整数型配列の読み込み
    subroutine get_data_array_int32_(unit, ret_array, ret_array_size)
        implicit none
        integer(4),intent(in) :: unit
        integer(4),allocatable,intent(out) :: ret_array(:)
        integer(4),intent(in) :: ret_array_size
        integer(4) ibyte, iretn, irecn
        integer(4) irec, L, ios
        integer(4) subrecn

        read(unit) ibyte, iretn, irecn

        subrecn = irecn - 1

        if(.not. allocated(ret_array)) allocate(ret_array(ret_array_size))

        if(subrecn == 0) then
            read(unit,iostat=ios) (ret_array(L), L=1, irecn*iretn)
            return
        else
            do irec = 1, subrecn
                read(unit,iostat=ios) (ret_array(L), L=1+(irec-1)*iretn, irec*iretn)
                if ( ios /= 0 ) then
                    ! print*,'iostat: ', ios, 'at L=', L
                    exit
                end if
            end do
        
            read(unit) (ret_array(L), L=1+(irecn-1)*iretn, ret_array_size)
        endif

    end subroutine

    !> 倍精度実数型配列の読み込み
    subroutine get_data_array_float64_(unit, ret_array, ret_array_size)
        implicit none
        integer(4),intent(in) :: unit
        real(8),allocatable,intent(out) :: ret_array(:)
        integer(4),intent(in) :: ret_array_size
        integer(4) ibyte, iretn, irecn
        integer(4) irec, L, ios
        integer(4) subrecn

        read(unit) ibyte, iretn, irecn

        subrecn = irecn - 1
        
        if(.not. allocated(ret_array)) allocate(ret_array(ret_array_size))

        if(subrecn == 0) then
            read(unit,iostat=ios) (ret_array(L), L=1, irecn*iretn)
            return
        else
            do irec = 1, subrecn
                read(unit,iostat=ios) (ret_array(L), L=1+(irec-1)*iretn, irec*iretn)
                if ( ios /= 0 ) then
                    ! print*,'iostat: ', ios, 'at L=', L
                    exit
                end if
            end do
        
            read(unit) (ret_array(L), L=1+(irecn-1)*iretn, ret_array_size)
        endif

    end subroutine

    !> データを読み飛ばす処理
    subroutine ignore_data_(unit)
        implicit none
        integer(4),intent(in) :: unit
        integer(4) ibyte, iretn, irecn
        integer(4) irec

        read(unit) ibyte, iretn, irecn
        do irec = 1, irecn
            read(unit)
        end do

    end subroutine

    !> 文字列（バイト数指定）
    subroutine get_data_char_(unit, byte, ret_char)
        implicit none
        integer(4), intent(in) :: unit
        integer(4), intent(in) :: byte
        character(:),allocatable,intent(inout) :: ret_char

        if(.not. allocated(ret_char)) allocate(character(byte) :: ret_char)
        read(unit) 
        read(unit) ret_char

    end subroutine

end module SCT_file_reader_m