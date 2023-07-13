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
        integer, allocatable :: adjacentCellIDs(:)
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
        integer,allocatable, private :: face2cells(:,:) 
        integer,allocatable, private :: cell2faces(:,:) 

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

        contains

        procedure, public :: read_SCF_file
        procedure, public :: get_fph_element_count
        procedure, public :: get_fph_vertex_count
        procedure, public :: get_fph_2d_array_of_point_coords
        procedure, public :: get_face2vertices
        procedure, public :: get_face2cells
        procedure, public :: get_fph_boundFaceIDs
        procedure, public :: output_fph_boundFace
        procedure, public :: get_fph_adjacentCellIDs
        procedure, public :: output_fph_adjacentCell

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
        class(scf_grid_t), intent(in) :: this

        get_fph_element_count = this%NELEM

    end function

    integer function get_fph_vertex_count(this)
        class(scf_grid_t), intent(in) :: this

        get_fph_vertex_count = this%NODES

    end function

    subroutine get_fph_2d_array_of_point_coords(this, points)
        !! 節点座標を2次元配列で出力する. 
        implicit none
        class(scf_grid_t),intent(in) :: this
        real(4), allocatable, intent(inout) :: points(:,:)

        allocate(points(3,this%NODES))

        points(1,:) = this%CAN_X(:)
        points(2,:) = this%CAN_Y(:)
        points(3,:) = this%CAN_Z(:)
        
    end subroutine

    subroutine get_face2vertices(this) 
        class(scf_grid_t), intent(inout) :: this
        integer ::  jj, kk, cnt

        allocate(this%face2vertices(this%NFACE))

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
        end do
            
    end subroutine

    subroutine get_face2cells(this) 
        class(scf_grid_t), intent(inout) :: this 
        integer :: jj

        allocate(this%face2cells(2,this%NFACE))

        ! セル番号を0番から1番スタートにする
        do jj = 1,this%NFACE
            this%face2cells(1,jj) = this%IE1(jj) + 1
            this%face2cells(2,jj) = this%IE2(jj) + 1   !0のときは存在しない 
        end do 

    end subroutine

    subroutine get_fph_boundFaceIDs(this)
        class(scf_grid_t), intent(inout) :: this
        integer jj
        logical first_flag

        first_flag = .true.

        ! セル番号0を有する面は境界面(外部表面)
        do jj = 1,this%NFACE
            if(this%face2cells(1,jj) == 0) this%boundFaceIDs = append2list_int(this%boundFaceIDs,jj,first_flag)
            if(this%face2cells(2,jj) == 0) this%boundFaceIDs = append2list_int(this%boundFaceIDs,jj,first_flag)
        end do

        this%num_boundFace = size(this%boundFaceIDs)

    end subroutine

    subroutine output_fph_boundFace(this, dir)
        class(scf_grid_t), intent(inout) :: this
        character(*), intent(in) :: dir
        integer JB, n_unit

        print*, 'OUTPUT:', dir//"boundary.txt"
        open(newunit = n_unit, file = dir//"boundary.txt" , status = 'replace')
            write(n_unit,*) this%num_boundFace
            do JB = 1, this%num_boundFace
                write(n_unit,'(*(g0:," "))') this%face2vertices(this%boundFaceIDs(JB))%vertexIDs
            end do
        close(n_unit)

    end subroutine

    subroutine get_fph_adjacentCellIDs(this)
        use terminalControler_m
        class(scf_grid_t), intent(inout) :: this
        integer ii, jj
        logical first_flag

        first_flag = .true.

        allocate(this%mainCell(this%NELEM))

        ! 計算コスト大,要改善
        print*, "Now solve adjacent cells ..."
        call set_formatTC('("Now solve fph adjacent cells ... [ #cellID : ",i6," / ",i6," ]")')
        do ii = 1, this%NELEM
            call print_progress([ii, this%NELEM])

            do jj = 1,this%NFACE
                if(this%face2cells(1,jj) == 0 .or. this%face2cells(2,jj) == 0) cycle
                if(this%face2cells(1,jj) == ii) &
                    this%mainCell(ii)%adjacentCellIDs &
                    = append2list_int(this%mainCell(ii)%adjacentCellIDs,this%face2cells(2,jj),first_flag)
                if(this%face2cells(2,jj) == ii) &
                    this%mainCell(ii)%adjacentCellIDs &
                    = append2list_int(this%mainCell(ii)%adjacentCellIDs,this%face2cells(1,jj),first_flag)
            end do
            
            first_flag = .true.
        end do

    end subroutine

    subroutine output_fph_adjacentCell(this,dir)
        class(scf_grid_t), intent(inout) :: this
        character(*), intent(in) :: dir
        integer n_unit, ii

        print*, 'OUTPUT:', dir//"adjacency.txt" 
        open(newunit = n_unit, file = dir//"adjacency.txt", status='replace')
            write(n_unit,'(*(g0:," "))') this%NELEM
            write(n_unit,'(*(g0:," "))') 100 !任意の数で大丈夫そう
            do ii = 1, this%NELEM
                write(n_unit, '(*(g0:," "))') size(this%mainCell(ii)%adjacentCellIDs), this%mainCell(ii)%adjacentCellIDs     
            end do
        close(n_unit)

    end subroutine

    function append2list_int(list, element, first_flag) result(after_list)
        integer, intent(in) :: list(:)
        integer, intent(in) :: element
        logical, intent(inout) :: first_flag
        integer, allocatable :: after_list(:)
        integer n

        if(first_flag) then
            allocate(after_list(1))
            after_list(1) = element
        else
            n = size(list)
            allocate(after_list(n+1))
            after_list(:n) = list(:n)
            after_list(n+1) = element
        end if

        first_flag = .false.

    end function

end module
