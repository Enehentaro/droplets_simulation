module SCF_file_reader_m
!================================================================================================================
!   version : 1.0.1     (2022/11/22)
!   author : Miyoshi 
!summary : 
!   SCFlow　出力ファイルを読み込み，データを取り出す
!   このプログラム単体で独立して扱えるようにする
!================================================================================================================

    implicit none

    real(8),parameter :: MissingValueSize = 1.0d20  !SCFlowでの欠測値の大きさ

    character(10),parameter :: EC_Scalar_HeadName = 'EC_Scalar:'
    character(10),parameter :: EC_Vector_HeadName = 'EC_Vector:'
    character(10),parameter :: FC_Scalar_HeadName = 'FC_Scalar:'
    character(10),parameter :: FC_Vector_HeadName = 'FC_Vector:'
    character(32),parameter :: OverlapEndLabel    = 'OverlapEnd' 
    
    
    integer,parameter :: EC_Scalar_DataSize = 50    !FPH内のセル重心の持つスカラーデータの個数上限（決め打ち）
    integer,parameter :: EC_Vector_DataSize = 50    !FPH内のセル重心の持つベクトルデータの個数上限（決め打ち）
    integer,parameter :: FC_Scalar_DataSize = 50    !FPH内の面重心の持つスカラーデータの個数上限（決め打ち）
    integer,parameter :: FC_Vector_DataSize = 50    !FPH内の面重心の持つベクトルデータの個数上限（決め打ち）

    type :: EC_Scalar_t 
        character(:),allocatable :: name 
        character(:),allocatable :: abbreviated_name 
        integer(4) ndata
        real(4),allocatable :: data(:) 
    end type

    type :: EC_Vector_t 
        character(:),allocatable :: name 
        character(:),allocatable :: abbreviated_name 
        integer(4) ndata  
        real(4),allocatable :: x(:),y(:),z(:)
    end type

    type :: FC_Scalar_t 
        character(:),allocatable :: name 
        character(:),allocatable :: abbreviated_name 
        integer(4) ndata
        integer(4),allocatable :: face_num(:) ,face_flag(:)
        real(4),allocatable :: data(:) 
    end type

    type :: FC_Vector_t 
        character(:),allocatable :: name 
        character(:),allocatable :: abbreviated_name 
        integer(4) ndata  
        integer(4),allocatable :: face_num(:),face_flag(:) 
        real(4),allocatable :: x(:),y(:),z(:)
    end type

    type :: data_name_list_t 
        character(:),allocatable :: name 
        character(:),allocatable :: abbreviated_name 
    end type

    type :: content_t
        integer, allocatable :: vertexIDs(:)
        integer, allocatable :: faceIDs(:)
            !! dummy (-99) を含む配列
        integer, allocatable :: boundFaceID(:)
            !! dummy (-99) を含む配列
        integer, allocatable :: adjacentCellIDs(:)
            !! dummy (-99) を含む配列
        real(4) center(3)
        real(4) coordinate(3)
    end type

    type :: scf_grid_t 
        
        real(4), private :: TIME = 0 
            !!時間　
        integer(4), private :: NCYC = 0 
            !!サイクル数

        integer(4), private :: NODES = 0                               
            !!節点数
        real(4), allocatable, private :: CAN_X(:),CAN_Y(:),CAN_Z(:)   
            !!節点座標
        integer(4), private :: NFACE = 0 
            !!要素界面数
        integer(4), private :: NELEM = 0
            !!要素数
        real(4),allocatable, private :: CCE_X(:),CCE_Y(:),CCE_Z(:)
            !!要素中心座標
        integer(4),allocatable, private :: IE1(:),IE2(:) 
            !!面の裏表の要素番号
        integer(4),allocatable, private :: NDNUM(:) 
            !!界面を構成する節点数
        integer(4), private :: NDTOT
            !!全界面を構成する節点数
        integer(4),allocatable, private :: IDNO(:) 
            !!界面を構成する節点番号

        type(content_t), allocatable, private :: face2vertices(:)
        type(content_t), allocatable, private :: mainCell(:)
        type(content_t), allocatable, private :: cell2faces(:)
        type(content_t), allocatable :: node(:)
        type(content_t), allocatable :: face(:)
        integer,allocatable, private :: face2cells(:,:) 

        integer, private :: EC_scalar_data_count = 0 
        integer, private :: EC_vector_data_count = 0 
        integer, private :: FC_scalar_data_count = 0 
        integer, private :: FC_vector_data_count = 0
        type(EC_Scalar_t),allocatable, private :: EC_Scalars(:) 
        type(EC_Vector_t),allocatable, private :: EC_Vectors(:) 
        type(FC_Scalar_t),allocatable, private :: FC_Scalars(:) 
        type(FC_Vector_t),allocatable, private :: FC_Vectors(:) 

        integer, private :: num_boundFace
        integer, allocatable, private :: boundFaceIDs(:)
            !! dummy (-99) を含む配列
        integer, allocatable, private :: num_face2vertex(:)
        integer, allocatable, private :: offsets(:)
        integer, private :: num_obj_vert
        integer, allocatable :: obj_pair_vertID(:)

        contains

        procedure, public :: read_SCF_file
        procedure, public :: get_fph_element_count
        procedure, public :: get_fph_vertex_count
        procedure, public :: get_fph_face_count
        procedure, public :: set_node_coords
        procedure, public :: get_fph_2d_array_of_point_coords
        procedure, public :: get_fph_2d_array_of_cell_coords
        procedure, public :: get_face2vertices
        procedure, public :: get_face2cells
        procedure, public :: set_cell2faces
        procedure, public :: get_cell2faces
        procedure, public :: get_fph_bound_faceIDs
        procedure, public :: get_fph_face_center
        procedure, public :: get_fph_bound_face_center
        procedure, public :: output_fph_cell2face
        procedure, public :: read_cell2face
        procedure, public :: output_fph_bound_face
        procedure, public :: drop_vertex_for_obj
        procedure, public :: output_OBJ
        procedure, public :: get_cell2bound_face
        procedure, public :: get_fph_adjacentCellIDs
        procedure, public :: output_fph_adjacentCell
        procedure, public :: search_fph_vector_data

    end type
    
    contains 

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

    subroutine read_SCF_file(this, filename)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        character(*),intent(in) :: filename
        integer unit

        !割り付けされている物があれば解放する.
        call destructor(this)

        if(.not.open_binary_sequential_(unit,filename))then 
            print*,'cannot open'
            stop 
        end if

        call read_FPH_Header_data(unit, this%NCYC, this%TIME)  
        call read_FPH_Main_data(unit, this%CAN_X, this%CAN_Y, this%CAN_Z, this%CCE_X, this%CCE_Y, this%CCE_Z, &
                                this%EC_Scalars, this%EC_Vectors, this%FC_Scalars, this%FC_Vectors, &
                                this%EC_scalar_data_count, this%EC_vector_data_count, &
                                this%FC_scalar_data_count, this%FC_vector_data_count, &
                                this%NODES, this%NFACE, this%NELEM, this%NDTOT, this%IE1, this%IE2, this%NDNUM, this%IDNO)

    end subroutine

    subroutine read_FPH_Header_data(unit,NCYC,TIME) 
        implicit none 
        integer,intent(in) :: unit
        integer(4), intent(inout) :: NCYC
        real(4), intent(inout) :: TIME
        character(8) header_name       
        integer(4) version_num, NNAMS , n     
        character(32) title_text 

        read(unit) header_name
        ! print*, 'header_name: ',header_name 
        read(unit) version_num
        ! print*, 'version_num :',version_num

        do 
            read(unit) title_text
            ! print*,title_text
            if(trim(title_text) == 'HeaderDataEnd') then 
                exit 
            else if (trim(title_text) == 'Cycle') then 
                call get_data_int32_(unit,NCYC) 
                call get_data_float64_(unit,TIME)
                read(unit)
                cycle 
            else if (trim(title_text) == 'Unused')then
                read(unit) 
                read(unit) 
            elseif(trim(title_text)=='Unit:$TEMP') then
                read(unit)
                read(unit)
                read(unit)
                read(unit)
                read(unit) NNAMS    !NNAMS : 温度とみなす変数の数
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

    subroutine read_FPH_Main_data(unit, CAN_X, CAN_Y, CAN_Z, CCE_X, CCE_Y, CCE_Z,&
                                    EC_Scalars, EC_Vectors, FC_Scalars, FC_Vectors ,&
                                    EC_Scalar_cnt, EC_Vector_cnt , FC_Scalar_cnt , FC_Vector_cnt,&
                                    NODES, NFACE, NELEM, NDTOT, IE1, IE2, NDNUM, IDNO)  

        implicit none 
        integer,intent(in) :: unit
        real(4), allocatable, intent(inout) :: CAN_X(:),CAN_Y(:),CAN_Z(:),CCE_X(:),CCE_Y(:),CCE_Z(:) 
        type(EC_Scalar_t), allocatable, intent(inout) :: EC_Scalars(:)
        type(EC_Vector_t), allocatable, intent(inout) :: EC_Vectors(:)
        type(FC_Scalar_t), allocatable, intent(inout) :: FC_Scalars(:)
        type(FC_Vector_t), allocatable, intent(inout) :: FC_Vectors(:)
        integer, intent(inout) :: EC_Scalar_cnt , EC_Vector_cnt , FC_Scalar_cnt , FC_Vector_cnt
        integer, intent(inout) :: NODES, NFACE, NELEM, NDTOT
        integer, allocatable, intent(inout) :: IE1(:), IE2(:), NDNUM(:), IDNO(:)
        integer :: NMAT, NMLEN, NPART, NPLEN, NREGN, NLEN, NN,  REV, NSIZE, LENG
        integer,allocatable :: NFA(:),FLG(:),ID(:),MAT(:),MAT_PART(:)
        integer i

        character(32) main_data_title 
        character(32) TITLE
        character,allocatable :: LMAT(:),LPART(:),LRGN(:),LRGN_S(:)
        character(:),allocatable ::  TEXT,STRXML


        

        read(unit) main_data_title  !'OverlapStart_n'<バイト数32>

        do 
            read(unit) TITLE
            ! print*,TITLE 

            if(trim(TITLE) == 'LS_Nodes' ) then 
                read(unit) !4,1,1
                read(unit) !LNX 
                call get_data_int32_(unit,NODES)                !接点数
                call get_data_array_float64_(unit,CAN_X,NODES)  !接点のX座標
                call get_data_array_float64_(unit,CAN_Y,NODES)  !接点のY座標
                call get_data_array_float64_(unit,CAN_Z,NODES)  !接点のZ座標
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'LS_Links') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NFACE)                !要素界面数　
                call get_data_array_int32_(unit,IE1,NFACE)      !裏面の要素番号
                call get_data_array_int32_(unit,IE2,NFACE)      !表面の要素番号（-1の場合、IE1の相手要素が存在しない）
                call get_data_array_int32_(unit,NDNUM,NFACE)    !面を構成する頂点数
                call get_data_int32_(unit,NDTOT)                !全界面を構成する接点数
                call get_data_array_int32_(unit,IDNO,NDTOT)     !界面を構成する接点番号（IE1からIE2へ向かう軸を右ねじ周りに構成）
                read(unit) !0,0,0 
            
            else if(trim(TITLE) =='LS_MaterialOfParts') then 
                read(unit) !4,1,1
                read(unit) !LNX 
                call get_data_int32_(unit,NMAT)                 !物性数
                call get_data_int32_(unit,NMLEN)                !LMATの文字数
                allocate(MAT(NMAT)) 
                allocate(LMAT(NMAT)) 
                do i = 1,NMAT 
                    read(unit) !4,1,1
                    read(unit) MAT(i)                           !物性番号
                    read(unit) !1,1,1 
                    read(unit) LMAT(i)                          !物性名 
                end do  
                call get_data_int32_(unit,NPART)                !部品数
                call get_data_int32_(unit,NPLEN)                !LPARTの文字数
                !allocate(LPART(NPART)) 
                allocate(MAT_PART(NPART))
                do i = 1,NPART 
                    read(unit) !1,1,1 
                    read(unit) LPART(i)                         !部品名
                    read(unit) !4,1,1 
                    read(unit) MAT_PART(i)                      !部品に割り当てられた物性番号
                end do 
                read (unit) !0,0,0 

            else if(trim(TITLE) == 'LS_CvolIdOfElements') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NELEM)                !要素数
                call get_data_array_int32_(unit,ID,NELEM)       !要素番号IEの閉空間ID
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'LS_SurfaceRegions') then 
                read(unit) !4,1,1 
                read(unit) !LNX  
                call get_data_int32_(unit,NREGN)                !領域データ数
                call get_data_int32_(unit,NLEN)                 !LRGNの文字数
                allocate(LRGN_S(NREGN))
                do i = 1,NREGN 
                    read(unit) !1,1,1 
                    read(unit) LRGN_S(i)                        !面領域名
                    read(unit) !4,1,1 
                    read(unit) NN                               !面領域を構成する面数
                    call get_data_array_int32_(unit,NFA,NN)     !面領域を構成する面番号
                    call get_data_array_int32_(unit,FLG,NN)     !面の裏表を指定するフラグ（1:IE1→IE2,2:IE2→IE1,3:IE1↔IE2）
                end do 
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'LS_VolumeRegions') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NREGN)                !領域データ数
                call get_data_int32_(unit,NLEN)                 !LRGNの文字数
                allocate(LRGN(NREGN)) 
                do i = 1,NREGN 
                    read(unit) !1,1,1 
                    read(unit) LRGN(i)                          !領域名
                    read(unit) !4,1,1 
                    read(unit) NN                               !体積領域を構成する閉空間数
                    call get_data_array_int32_(unit,ID,NN)      !体積領域を構成する閉空間ID
                end do 
                read(unit) !0,0,0

            else if(trim(TITLE) == 'LS_Parts') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NPART)                !部品データ数
                call get_data_int32_(unit,NLEN)                 !LPRTの文字数
                allocate(LPART(NPART)) 
                do i = 1,NPART 
                    read(unit) !1,1,1 
                    read(unit) LPART(i)                         !部品名
                    read(unit) !4,4,4 
                    read(unit) NN                               !部品を構成する閉空間数
                    call get_data_array_int32_(unit,ID,NN)      !部品を構成する閉空間ID
                end do 
                read(unit) !0,0,0  

            else if(trim(TITLE) == 'LS_Assemblies') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NLEN)                 !STRXMLのバイト数
                call get_data_char_(unit,NLEN,STRXML)           !部品とアセンブリの関する情報（xml形式）
                ! print*,STRXML
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'LS_PartialInformation') then 
                read(unit) !4,1,1 
                read(unit) !LNX  
                call get_data_int32_(unit,REV)                  !リビジョン番号1を指定
                call get_data_int32_(unit,NSIZE)                !TEXTのバイト単位でのサイズ
                call get_data_char_(unit,NSIZE,TEXT)            !部分図化ファイルの位置情報（無視しても良い?）
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'Element_Center') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,NELEM)                !フィールド変数のデータ数[=要素数]
                call get_data_array_float64_(unit,CCE_X,NELEM)  !要素中心のX座標
                call get_data_array_float64_(unit,CCE_Y,NELEM)  !要素中心のY座標
                call get_data_array_float64_(unit,CCE_Z,NELEM)  !要素中心のZ座標
                read(unit) !0,0,0 

            else if(trim(TITLE) == 'LS_SPHFile') then 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_int32_(unit,LENG)                 !TEXTの文字数
                call get_data_char_(unit,LENG,TEXT)             !SPHファイルの中身のテキスト（いらん?）
                ! print*,TEXT
                read(unit) !0,0,0 

            else if(TITLE(1:10) == 'EC_Scalar:') then 
                !EC_Scalar:以下はループ外で処理するので、1行戻ってループを抜ける
                backspace(unit)
                exit
                
            else if(TITLE(1:10) == 'EC_Vector:') then 
                !EC_Scalar:以下はループ外で処理するので、1行戻ってループを抜ける
                backspace(unit)
                exit

            else
                print*,'(error) readFPH_main_data :: UnKnown title data was detected.' , trim(TITLE)
                return 
            end if 
        end do 

        if(.not.allocated(EC_Scalars)) allocate(EC_Scalars(EC_Scalar_DataSize)) 
        if(.not.allocated(EC_Vectors)) allocate(EC_Vectors(EC_Vector_DataSize)) 
        if(.not.allocated(FC_Scalars)) allocate(FC_Scalars(FC_Scalar_DataSize)) 
        if(.not.allocated(FC_Vectors)) allocate(FC_Vectors(FC_Vector_DataSize)) 
        
        datas_reader: block 
            logical is_end 
            integer EC_s_cnt , EC_v_cnt , FC_s_cnt , FC_v_cnt , n

            is_end = .false. 
            EC_v_cnt = 0 
            do while(.not.is_end) 
                EC_v_cnt = EC_v_cnt + 1 
                call EC_Vector_reader(EC_Vectors(EC_v_cnt),unit,is_end) 
            end do 
            EC_Vector_cnt = EC_v_cnt - 1 

            is_end = .false. 
            EC_s_cnt = 0 
            do while(.not.is_end)   !is_endが.ture.になるまでdoループが回る。
                EC_s_cnt = EC_s_cnt + 1 
                call EC_Scalar_reader(EC_Scalars(EC_s_cnt),unit,is_end)
            end do
            EC_Scalar_cnt = EC_s_cnt - 1 

            is_end = .false. 
            FC_v_cnt = 0 
            do while(.not.is_end) 
                FC_v_cnt = FC_v_cnt + 1 
                call FC_Vector_reader(FC_Vectors(FC_v_cnt),unit,is_end) 
            end do 
            FC_Vector_cnt = FC_v_cnt - 1 

            is_end = .false. 
            FC_s_cnt = 0 
            do while(.not.is_end)  
                FC_s_cnt = FC_s_cnt + 1 
                call FC_Scalar_reader(FC_Scalars(FC_s_cnt),unit,is_end)
            end do
            FC_Scalar_cnt = FC_s_cnt - 1 

            
            do n = EC_s_cnt,size(EC_Scalars) 
                EC_Scalars(n)%name = 'NONE' 
            end do 
            
            do n = EC_v_cnt,size(EC_Vectors) 
                EC_Vectors(n)%name = 'NONE' 
            end do 
            do n = FC_s_cnt,size(FC_Scalars) 
                FC_Scalars(n)%name = 'NONE' 
            end do 
            
            do n = FC_v_cnt,size(FC_Vectors) 
                FC_Vectors(n)%name = 'NONE' 
            end do 
        end block datas_reader 

    end subroutine


    subroutine EC_Scalar_reader(scalar,unit,is_end)
        implicit none 
        type(EC_Scalar_t) :: scalar
        integer(4),intent(in) :: unit 
        logical,intent(out) :: is_end 
        character(32) title 

        is_end = .false. 
        read(unit) title
        ! print*,title
        select case(title(1:10)) 
            case(EC_Scalar_HeadName)
                scalar%abbreviated_name = title(11:14) 
                read(unit) !4,1,1 
                read(unit) !LNX  
                call get_data_char_(unit,32,scalar%name) 
                call get_data_int32_(unit,scalar%ndata) 
                call get_data_array_float64_(unit,scalar%data,scalar%ndata) 
                read(unit) !0,0,0 
            case(EC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(FC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(FC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(OverlapEndLabel) 
                is_end = .true. 
                return
            case default 
                print*,'something is wrong.',trim(title(1:10)) 
        end select

    end subroutine


    subroutine EC_Vector_reader(vector,unit,is_end)
        implicit none 
        type(EC_Vector_t) :: vector 
        integer(4),intent(in) :: unit 
        logical,intent(out) :: is_end 
        character(32) title 

        is_end = .false. 
        read(unit) title 
        ! print*,title
        select case(title(1:10))
            case(EC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(EC_Vector_HeadName)
                vector%abbreviated_name = title(11:14) 
                read(unit) !4,1,1 
                read(unit) !LNX
                call ignore_data_(unit) !LVCT 位置ベクトルなら1・そうでないなら0　必要ないと思うので無視
                call get_data_char_(unit,32,vector%name) 
                call get_data_int32_(unit,vector%ndata) 
                call get_data_array_float64_(unit,vector%x,vector%ndata)
                call get_data_array_float64_(unit,vector%y,vector%ndata)
                call get_data_array_float64_(unit,vector%z,vector%ndata)
                read(unit) !0,0,0 
            case(FC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(FC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(OverlapEndLabel)
                is_end = .true. 
                return 
            case default 
                print*,'something is wrong.',trim(title(1:10)) 
        end select

    end subroutine

    subroutine FC_Scalar_reader(scalar,unit,is_end) 
        implicit none 
        type(FC_Scalar_t) :: scalar 
        integer(4),intent(in) :: unit 
        logical,intent(out):: is_end 
        character(32) title 

        is_end = .false. 
        read(unit) title 
        ! print*,title
        select case(title(1:10))
            case(EC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return  
            case(EC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return  
            case(FC_Scalar_HeadName) 
                scalar%abbreviated_name = title(11:14) 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call get_data_char_(unit,32,scalar%name) 
                call get_data_int32_(unit,scalar%ndata) 
                call get_data_array_int32_(unit,scalar%face_num,scalar%ndata)
                call get_data_array_int32_(unit,scalar%face_flag,scalar%ndata) 
                call get_data_array_float64_(unit,scalar%data,scalar%ndata) 
                read(unit) !0,0,0 
            case(FC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return 
            case(OverlapEndLabel)
                is_end = .true. 
                return 
            case default 
                print*,'something is wrong.',trim(title(1:10)) 
        end select 
    end subroutine 

    
    subroutine FC_Vector_reader(vector,unit,is_end) 
        implicit none 
        type(FC_Vector_t) :: vector 
        integer(4),intent(in) :: unit 
        logical,intent(out) :: is_end 
        character(32) title 

        is_end = .false. 
        read(unit) title 
        ! print*,title
        select case(title(1:10)) 
            case(EC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return  
            case(EC_Vector_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return  
            case(FC_Scalar_HeadName) 
                backspace(unit) 
                is_end = .true. 
                return  
            case(FC_Vector_HeadName)
                vector%abbreviated_name = title(11:14) 
                read(unit) !4,1,1 
                read(unit) !LNX 
                call ignore_data_(unit) !LVCT 位置ベクトルなら1・そうでないなら0　必要ないと思うので無視
                call get_data_char_(unit,32,vector%name) 
                call get_data_int32_(unit,vector%ndata) 
                call get_data_array_int32_(unit,vector%face_num,vector%ndata) 
                call get_data_array_int32_(unit,vector%face_flag,vector%ndata) 
                call get_data_array_float64_(unit,vector%x,vector%ndata)
                call get_data_array_float64_(unit,vector%y,vector%ndata)
                call get_data_array_float64_(unit,vector%z,vector%ndata)
                read(unit) !0,0,0 
            case(OverlapEndLabel)
                is_end = .true. 
                return 
            case default 
                print*,'something is wrong.',trim(title(1:10)) 
                error stop 
        end select 

    end subroutine

    subroutine get_data_int32_(unit,retval)
        implicit none
        integer(4),intent(in) :: unit
        integer(4),intent(inout) :: retval 
        integer(4) ibyte, iretn, irecn 

        read(unit) ibyte, iretn, irecn
        read(unit) retval 
        
    end subroutine

    subroutine get_data_float64_(unit,retval)
        implicit none
        integer(4),intent(in) :: unit
        real(4),intent(inout) :: retval 
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
        real(4),allocatable,intent(out) :: ret_array(:)
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
                    print*,'iostat: ', ios, 'at L=', L
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

    subroutine destructor(this)
        implicit none
        type(scf_grid_t),intent(inout) :: this
    
        if (allocated(this%CAN_X))      deallocate(this%CAN_X)
        if (allocated(this%CAN_Y))      deallocate(this%CAN_Y)
        if (allocated(this%CAN_Z))      deallocate(this%CAN_Z)
        if (allocated(this%CCE_X))      deallocate(this%CCE_X)
        if (allocated(this%CCE_Y))      deallocate(this%CCE_Y)
        if (allocated(this%CCE_Z))      deallocate(this%CCE_Z)
        if (allocated(this%EC_Scalars)) deallocate(this%EC_Scalars)
        if (allocated(this%EC_Vectors)) deallocate(this%EC_Vectors)
        if (allocated(this%FC_Scalars)) deallocate(this%FC_Scalars)
        if (allocated(this%EC_Vectors)) deallocate(this%EC_Vectors)
        if (allocated(this%face2vertices)) deallocate(this%face2vertices)
    
        this%NODES = 0
        this%NFACE = 0 
        this%NELEM = 0
        this%EC_scalar_data_count = 0 
        this%EC_vector_data_count = 0 
        this%FC_scalar_data_count = 0 
        this%FC_vector_data_count = 0 
    end subroutine

    integer function get_fph_element_count(this)
        implicit none
        class(scf_grid_t), intent(in) :: this

        get_fph_element_count = this%NELEM

    end function

    integer function get_fph_vertex_count(this)
        implicit none
        class(scf_grid_t), intent(in) :: this

        get_fph_vertex_count = this%NODES

    end function

    integer function get_fph_face_count(this)
        implicit none
        class(scf_grid_t), intent(in) :: this

        get_fph_face_count = this%NFACE

    end function

    subroutine get_fph_2d_array_of_point_coords(this, points)
        !! 節点座標を2次元配列で出力する. 
        implicit none
        class(scf_grid_t),intent(inout) :: this
        real(4), allocatable, intent(inout) :: points(:,:)

        call packing_vector_into_2Darray_(points, this%CAN_X, this%CAN_Y, this%CAN_Z)

    end subroutine

    subroutine get_fph_2d_array_of_cell_coords(this, cells)
        !! 要素中心座標を2次元配列で出力する. 
        implicit none
        class(scf_grid_t),intent(inout) :: this
        real(4), allocatable, intent(inout) :: cells(:,:)

        call packing_vector_into_2Darray_(cells, this%CCE_X, this%CCE_Y, this%CCE_Z)
        
    end subroutine

    subroutine set_node_coords(this)
        class(scf_grid_t),intent(inout) :: this
        integer kk

        allocate(this%node(this%NODES))

        do kk = 1, this%NODES
            this%node(kk)%coordinate(1) = this%CAN_X(kk)
            this%node(kk)%coordinate(2) = this%CAN_Y(kk)
            this%node(kk)%coordinate(3) = this%CAN_Z(kk)
        end do

    end subroutine

    subroutine get_face2vertices(this)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        integer ::  jj, kk, cnt

        allocate(this%face2vertices(this%NFACE))
        allocate(this%num_face2vertex(this%NFACE))

        do jj = 1, this%NFACE
            allocate(this%face2vertices(jj)%vertexIDs(this%NDNUM(jj)))
        end do

        cnt = 1
        ! 頂点番号を0番スタートから1番スタートにする
        do jj = 1, this%NFACE
            do kk = 1, this%NDNUM(jj)
                this%face2vertices(jj)%vertexIDs(kk) = this%IDNO(cnt) + 1
                cnt = cnt + 1
            end do
            this%num_face2vertex(JJ) = size(this%face2vertices(jj)%vertexIDs)        
        end do
            
    end subroutine

    subroutine get_face2cells(this)
        implicit none
        class(scf_grid_t), intent(inout) :: this 
        integer :: jj

        allocate(this%face2cells(2,this%NFACE))

        ! セル番号を0番から1番スタートにする
        do jj = 1,this%NFACE
            this%face2cells(1,jj) = this%IE1(jj) + 1
            this%face2cells(2,jj) = this%IE2(jj) + 1   !0のときは存在しない 
        end do 

    end subroutine

    subroutine get_fph_face_center(this,face_center)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        real(4), allocatable, intent(out) :: face_center(:,:)
        real(4) sum_vertices(3)
        integer jj, vertexID

        allocate(face_center(3,this%NFACE))
        allocate(this%face(this%NFACE))

        do jj = 1, this%NFACE

            sum_vertices(:) = 0.0
            do vertexID = 1, size(this%face2vertices(jj)%vertexIDs)
                sum_vertices = sum_vertices + &
                this%node(this%face2vertices(jj)%vertexIDs(vertexID))%coordinate
            end do

            this%face(jj)%center(:) = sum_vertices / size(this%face2vertices(jj)%vertexIDs)
            face_center(:,jj) = this%face(jj)%center(:)
            
        end do
        
    end subroutine

    subroutine get_fph_bound_faceIDs(this, num_boundFaces)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        integer, intent(out) :: num_boundFaces
        integer jj, kk, ll, alloc_max, dummyID

        ! 配列のサイズが未確定なのでダミー配列を用意
        alloc_max = this%NFACE
        allocate(this%boundFaceIDs(alloc_max), source = -99)

        ! セル番号0を有する面は境界面(外部表面)
        do jj = 1,this%NFACE
            if(any(this%face2cells(:,jj) == 0)) then
                do kk = 1, this%NFACE
                    if(this%boundFaceIDs(kk) == -99) then
                        dummyID = kk
                        exit
                    endif
                end do
                ! dummyID = findloc(this%boundFaceIDs, -99, dim = 1)
                call check_range_of_array(dummyID, alloc_max, "boundFaceIDs")
                this%boundFaceIDs(dummyID) = jj
            end if
        end do

        ! -99が初めて見つかった時の番地をdummyIDに格納
        do ll = 1, this%NFACE
            if(this%boundFaceIDs(ll) == -99) then
                dummyID = ll
                exit
            endif
        end do
        ! dummyID = findloc(this%boundFaceIDs, -99, dim = 1)

        this%num_boundFace = dummyID-1
        num_boundFaces = this%num_boundFace

    end subroutine

    subroutine get_fph_bound_face_center(this, bound_center)
        implicit none
        class(scf_grid_t), intent(in) :: this
        real(4), allocatable, intent(inout) :: bound_center(:,:)
        integer JB

        allocate(bound_center(3,this%num_boundFace))

        do JB = 1, this%num_boundFace
            bound_center(:,JB) = this%face(this%boundFaceIDs(JB))%center(:)
        end do

    end subroutine

    subroutine set_cell2faces(this)
        use terminalControler_m
        !$ use omp_lib
        implicit none
        class(scf_grid_t), intent(inout) :: this
        integer jj, cellID, faceID, contentID, alloc_max, dummyID

        allocate(this%cell2faces(this%NELEM))
        ! 配列のサイズが未確定なのでダミー配列を用意
        alloc_max = 100
        do cellID = 1, this%NELEM
            allocate(this%cell2faces(cellID)%faceIDs(alloc_max), source = -99)
        end do
        
        print*, "Now get cell2face ..."
        call set_formatTC('("completed ... [ #faceID : ",i8," / ",i8," ]")')

        !$omp parallel do
        do faceID = 1, this%NFACE
            call print_progress([faceID, this%NFACE])
            do contentID = 1,2
                cellID = this%face2cells(contentID,faceID)

                ! cellIDが0のときはスルー
                if(cellID == 0) cycle

                ! 左から数えて何番目に-99があるか探索
                do jj = 1, alloc_max
                    if(this%cell2faces(cellID)%faceIDs(jj) == -99) then
                        dummyID = jj
                        exit
                    endif
                end do
                ! dummyID = findloc(this%cell2faces(cellID)%faceIDs, -99, dim = 1)

                call check_range_of_array(dummyID, alloc_max, "cell2faces")

                ! -99をfaceIDに置き換え
                this%cell2faces(cellID)%faceIDs(dummyID) = faceID
            end do
        end do
        !$omp end parallel do

    end subroutine

    subroutine output_fph_cell2face(this, dir)
        implicit none
        class(scf_grid_t), intent(in) :: this
        character(*), intent(in) :: dir
        integer n_unit, cellID

        print*, 'OUTPUT:', dir//"cell2face.txt"
        open(newunit = n_unit, file = dir//"cell2face.txt", status = 'replace')
            write(n_unit, '(*(g0:," "))') size(this%cell2faces(1)%faceIDs)
            do cellID = 1, this%NELEM
                write(n_unit, '(*(g0:," "))') this%cell2faces(cellID)%faceIDs
            end do
        close(n_unit)        

    end subroutine

    subroutine read_cell2face(this, dir)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        character(*), intent(in) :: dir
        integer n_unit, num_cell2face, cellID

        allocate(this%cell2faces(this%NELEM))
        
        open(newunit = n_unit, file = dir//'cell2face.txt', status = 'old')
            read(n_unit,*) num_cell2face
            do cellID = 1, this%NELEM
                allocate(this%cell2faces(cellID)%faceIDs(num_cell2face))
                read(n_unit, *) this%cell2faces(cellID)%faceIDs(:)
            end do
        close(n_unit)

    end subroutine

    function get_cell2faces(this) result(cell2face)
        implicit none
        class(scf_grid_t), intent(in) :: this
        integer, allocatable :: cell2face(:,:)
        integer cellID

        allocate(cell2face(this%NELEM,size(this%cell2faces(1)%faceIDs)))
        do cellID = 1, this%NELEM
            cell2face(cellID,:) = this%cell2faces(cellID)%faceIDs(:)
        end do

    end function

    subroutine get_cell2bound_face(this)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        integer bFID_1, bFID_2, JB, boundFace2cellID, alloc_max, cellID, dummyID

        if(.not.allocated(this%mainCell)) allocate(this%mainCell(this%NELEM))
        ! 配列のサイズが未確定なのでダミー配列を用意
        alloc_max = 10
        do cellID = 1, this%NELEM
            allocate(this%mainCell(cellID)%boundFaceID(alloc_max), source = -99)
        end do

        ! 境界セルに境界面番号を割り当てる
        do JB = 1, this%num_boundFace

            if(this%face2cells(1,this%boundFaceIDs(JB)) == 0) then

                boundFace2cellID = this%face2cells(2,this%boundFaceIDs(JB))
                do bFID_1 = 1, alloc_max
                    if(this%mainCell(boundFace2cellID)%boundFaceID(bFID) == -99) then
                        dummyID = bFID
                        exit
                    endif
                end do
                ! dummyID = findloc(this%mainCell(boundFace2cellID)%boundFaceID, -99, dim = 1)
                call check_range_of_array(dummyID, alloc_max, "mainCell%boundFace")
                this%mainCell(boundFace2cellID)%boundFaceID(dummyID) = JB

            end if
            if(this%face2cells(2,this%boundFaceIDs(JB)) == 0) then
                
                boundFace2cellID = this%face2cells(1,this%boundFaceIDs(JB))
                do bFID_2 = 1, alloc_max
                    if(this%mainCell(boundFace2cellID)%boundFaceID(bFID_2) == -99) then
                        dummyID = bFID_2
                        exit
                    endif
                end do
                ! dummyID = findloc(this%mainCell(boundFace2cellID)%boundFaceID, -99, dim = 1)
                call check_range_of_array(dummyID, alloc_max, "mainCell%boundFace")
                this%mainCell(boundFace2cellID)%boundFaceID(dummyID) = JB

            end if
        end do

        ! 未割り当ての配列は内部セル
        do cellID = 1, this%NELEM
            if(all(this%mainCell(cellID)%boundFaceID == -99)) then
                this%mainCell(cellID)%boundFaceID(1) = 0
            end if
        end do     

    end subroutine

    subroutine get_fph_adjacentCellIDs(this)
        use terminalControler_m
        !$ use omp_lib
        implicit none
        class(scf_grid_t), intent(inout) :: this
        integer faceID, cellID, alloc_max, dummyID, mainCellID, adjacentCellID, adID_1, adID_2

        ! 配列のサイズが未確定なのでダミー配列を用意
        alloc_max = 100
        do cellID = 1, this%NELEM
            allocate(this%mainCell(cellID)%adjacentCellIDs(alloc_max), source = -99)
        end do

        ! 計算コスト大,要改善
        print*, "Now solve adjacent cells ..."
        call set_formatTC('("Completed ... [ #faceID : ",i8," / ",i8," ]")')

        !$omp parallel do
        do faceID = 1, this%NFACE
            call print_progress([faceID, this%NFACE])

            if(all(this%face2cells(:,faceID) /= 0)) then
                mainCellID = this%face2cells(1,faceID)
                adjacentCellID = this%face2cells(2,faceID)
                do adID_1 = 1, alloc_max
                    if(this%mainCell(mainCellID)%adjacentCellIDs(adID_1) == -99) then
                        dummyID = adID_1
                        exit
                    end if
                end do
                ! dummyID = findloc(this%mainCell(mainCellID)%adjacentCellIDs, -99, dim = 1)
                call check_range_of_array(dummyID, alloc_max, "mainCell%adjacentCellIDs")
                this%mainCell(mainCellID)%adjacentCellIDs(dummyID) = adjacentCellID

                mainCellID = this%face2cells(2,faceID)
                adjacentCellID = this%face2cells(1,faceID)
                do adID_2 = 1, alloc_max
                    if(this%mainCell(mainCellID)%adjacentCellIDs(adID_2) == -99) then
                        dummyID = adID_2
                        exit
                    end if
                end do
                ! dummyID = findloc(this%mainCell(mainCellID)%adjacentCellIDs, -99, dim = 1)
                call check_range_of_array(dummyID, alloc_max, "mainCell%adjacentCellIDs")
                this%mainCell(mainCellID)%adjacentCellIDs(dummyID) = adjacentCellID
            end if

        end do
        !$omp end parallel do

    end subroutine

    subroutine output_fph_bound_face(this, dir)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        character(*), intent(in) :: dir
        integer JB, n_unit

        print*, 'OUTPUT:', dir//"boundary.txt"
        open(newunit = n_unit, file = dir//"boundary.txt" , status = 'replace')
            write(n_unit,'(i0)') this%num_boundFace
            do JB = 1, this%num_boundFace
                write(n_unit,'(*(g0:," "))') size(this%face2vertices(this%boundFaceIDs(JB))%vertexIDs) ,&
                this%face2vertices(this%boundFaceIDs(JB))%vertexIDs
            end do
        close(n_unit)

    end subroutine

    subroutine output_fph_adjacentCell(this,dir)
        implicit none
        class(scf_grid_t), intent(inout) :: this
        character(*), intent(in) :: dir
        integer n_unit, cellID, dummyID, adID, alloc_max_ad, alloc_max_bf, bfID

        alloc_max_ad = 100
        alloc_max_bf = 10
        print*, 'OUTPUT:', dir//"adjacency.txt" 
        open(newunit = n_unit, file = dir//"adjacency.txt", status='replace')
            write(n_unit,'(*(g0:," "))') this%NELEM
            write(n_unit,'(*(g0:," "))') 100 !任意の数で大丈夫そう
            do cellID = 1, this%NELEM
                do adID = 1, alloc_max_ad
                    if(this%mainCell(cellID)%adjacentCellIDs(adID) == -99) then
                        dummyID = adID
                        exit
                    end if
                end do
                ! dummyID = findloc(this%mainCell(cellID)%adjacentCellIDs, -99, dim = 1)
                write(n_unit, '(*(g0:," "))') dummyID-1, this%mainCell(cellID)%adjacentCellIDs(1:dummyID-1)
            end do
            do cellID = 1, this%NELEM
                if(this%mainCell(cellID)%boundFaceID(1) == 0) then
                    write(n_unit, '(*(g0:," "))') 0
                else
                    do bfID = 1, alloc_max_bf
                        if(this%mainCell(cellID)%boundFaceID(bfID) == -99) then
                            dummyID = bfID
                            exit
                        end if
                    end do
                    ! dummyID = findloc(this%mainCell(cellID)%boundFaceID, -99, dim = 1)
                    write(n_unit, '(*(g0:," "))') dummyID-1, this%mainCell(cellID)%boundFaceID(1:dummyID-1)                
                end if
            end do
        close(n_unit)

    end subroutine

    
    subroutine drop_vertex_for_obj(this)
        class(scf_grid_t), intent(inout) :: this
        integer bd_faceID, elementID, bd_vertID, cnt, JB

        allocate(this%obj_pair_vertID(this%NODES), source = -1)

        cnt = 1
        do JB = 1, this%num_boundFace
            bd_faceID = this%boundFaceIDs(JB)
            do elementID = 1, this%num_face2vertex(bd_faceID)
                bd_vertID = this%face2vertices(bd_faceID)%vertexIDs(elementID)
                if(this%obj_pair_vertID(bd_vertID) == -1) then
                    this%obj_pair_vertID(bd_vertID) = cnt
                    cnt = cnt + 1
                end if
            end do
        end do

        this%num_obj_vert = cnt

    end subroutine

    subroutine output_OBJ(this, dir)
        class(scf_grid_t), intent(in) :: this

        character(*), intent(in) :: dir
        integer n_unit, bd_faceID, elementID, bd_vertID, JB!, nodeID
        logical, allocatable :: is_written(:)

        allocate(is_written(this%NODES), source = .false.)

        print*, 'OUTPUT:', dir//"shape.obj"
        open(newunit = n_unit, file = dir//"shape.obj", status = "replace")
            write(n_unit, "(a)") "# Vertices"
            ! do nodeID = 1, this%NODES
            !     write(n_unit, '(*(g0:," "))') "v", this%node(nodeID)%coordinate
            ! end do
            do JB = 1, this%num_boundFace
                bd_faceID = this%boundFaceIDs(JB)
                do elementID = 1, this%num_face2vertex(bd_faceID)
                    bd_vertID = this%face2vertices(bd_faceID)%vertexIDs(elementID)
                    if(.not.is_written(bd_vertID)) then
                        write(n_unit, '(*(g0:," "))') "v", this%node(bd_vertID)%coordinate
                        is_written(bd_vertID) = .true.
                    end if
                end do
            end do
            write(n_unit, "(a)") "# Faces"
            do JB = 1, this%num_boundFace
                bd_faceID = this%boundFaceIDs(JB)
                write(n_unit, '(a)', advance = "no") "f "
                do elementID = 1, this%num_face2vertex(bd_faceID)
                    bd_vertID = this%face2vertices(bd_faceID)%vertexIDs(elementID)
                    write(n_unit, '(i0,1x)', advance = "no") this%obj_pair_vertID(bd_vertID)
                end do
                write(n_unit,"()")
            end do
        close(n_unit)

    end subroutine

    subroutine search_fph_vector_data(this, key, vector)
        implicit none
        class(scf_grid_t), intent(in) :: this
        character(*),intent(in) :: key
        real(4), allocatable, intent(inout) :: vector(:,:)
        integer i

        do i = 1, size(this%EC_Vectors)
            if ( trim(this%EC_Vectors(i)%abbreviated_name) == trim(key) ) then
                call packing_vector_into_2Darray_(vector, &
                this%EC_Vectors(i)%x, this%EC_Vectors(i)%y, this%EC_Vectors(i)%z)
                return
            end if
        end do

    end subroutine

    subroutine packing_vector_into_2Darray_(array, x, y, z)
        !! 実数のベクトル配列を2次元配列に詰め直す
        implicit none
        real(4), allocatable, intent(out) :: array(:,:)
        real(4), intent(in) :: x(:), y(:), z(:)

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

    subroutine check_range_of_array(allocID, alloc_max, array_name)
        integer, intent(in) :: allocID, alloc_max
        character(*), intent(in) :: array_name

        if(allocID > alloc_max) then
            print*, "allocate beyond upper bound of "//array_name
            stop
        end if

    end subroutine

end module
