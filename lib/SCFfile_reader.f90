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

    type :: scf_grid_t 
        
        real(4) :: TIME = 0 
            !!時間　
        integer(4) :: NCYC = 0 
            !!サイクル数

        integer(4) :: NODES = 0                               
            !!節点数
        real(4),allocatable :: CAN_X(:),CAN_Y(:),CAN_Z(:)   
            !!節点座標
        integer(4) :: NFACE = 0 
            !!要素界面数
        integer(4) :: NELEM = 0
            !!要素数
        real(4),allocatable :: CCE_X(:),CCE_Y(:),CCE_Z(:)
            !!要素中心座標
        integer(4),allocatable :: IE1(:),IE2(:) 
            !!面の裏表の要素番号
        integer(4),allocatable :: NDNUM(:) 
            !!界面を構成する節点数
        integer(4) :: NDTOT
            !!全界面を構成する節点数
        integer(4),allocatable :: IDNO(:) 
            !!界面を構成する節点番号


        integer,allocatable :: face2vertices(:,:)
        integer,allocatable :: face2cells(:,:) 
        integer,allocatable :: cell2faces(:,:) 

        integer :: EC_scalar_data_count = 0 
        integer :: EC_vector_data_count = 0 
        integer :: FC_scalar_data_count = 0 
        integer :: FC_vector_data_count = 0
        type(EC_Scalar_t),allocatable :: EC_Scalars(:) 
        type(EC_Vector_t),allocatable :: EC_Vectors(:) 
        type(FC_Scalar_t),allocatable :: FC_Scalars(:) 
        type(FC_Vector_t),allocatable :: FC_Vectors(:) 

        character(len=20) case_name 

        contains

        procedure :: read_SCF_file

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
                                this%NODES, this%NFACE, this%NELEM ,this%EC_scalar_data_count, this%EC_vector_data_count, &
                                this%FC_scalar_data_count, this%FC_vector_data_count, &
                                this%IE1, this%IE2, this%NDNUM, this%NDTOT, this%IDNO)

    end subroutine

    subroutine read_FPH_Header_data(unit,NCYC,TIME) 
        implicit none 
        integer,intent(in) :: unit 
        character(8) header_name       
        integer(4) version_num , NCYC, NNAMS , n     
        real(4) TIME 
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

    subroutine read_FPH_Main_data(unit, CAN_X, CAN_Y, CAN_Z, CCE_X, CCE_Y, CCE_Z,EC_Scalars, EC_Vectors, FC_Scalars, FC_Vectors ,&
                                    NODES, NFACE, NELEM, EC_Scalar_cnt, EC_Vector_cnt , FC_Scalar_cnt , FC_Vector_cnt, & 
                                    IE1, IE2, NDNUM, NDTOT, IDNO)  

        implicit none 
        integer,intent(in) :: unit 
        integer :: EC_Scalar_cnt , EC_Vector_cnt , FC_Scalar_cnt , FC_Vector_cnt, i
        integer :: NODES, NFACE, NDTOT, NMAT, NMLEN, NPART, NPLEN, NREGN, NLEN, NN, NELEM, REV, NSIZE, LENG
        integer,allocatable :: IE1(:),IE2(:),NDNUM(:),NFA(:),FLG(:),ID(:),MAT(:),MAT_PART(:),IDNO(:)
        real(4),allocatable :: CAN_X(:),CAN_Y(:),CAN_Z(:),CCE_X(:),CCE_Y(:),CCE_Z(:) 
        character(32) main_data_title 
        character(32) TITLE
        character,allocatable :: LMAT(:),LPART(:),LRGN(:),LRGN_S(:)
        character(:),allocatable ::  TEXT,STRXML
        type(EC_Scalar_t),allocatable :: EC_Scalars(:)
        type(EC_Vector_t),allocatable :: EC_Vectors(:)
        type(FC_Scalar_t),allocatable :: FC_Scalars(:)
        type(FC_Vector_t),allocatable :: FC_Vectors(:)

        

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

    subroutine make_face2vertices(this) 
        type(scf_grid_t) :: this 
        integer ::  jj, jjmx, kk,kkmx ,max_vertices, cnt_vertices, cnt 

        jjmx = this%NFACE 
        kkmx = this%NODES 
        max_vertices = 0 
        cnt_vertices = 0 
        do jj = 1,jjmx 
            cnt_vertices = this%NDNUM(jj) 
            if(cnt_vertices > max_vertices)then 
                max_vertices = cnt_vertices
            end if 
        end do
        ! print*,'max_vertices in a face is',max_vertices

        allocate(this%face2vertices(max_vertices,jjmx)) 
        
        cnt = 0 
        do jj = 1,jjmx 
            cnt_vertices = this%NDNUM(jj) 
            do kk = 1,cnt_vertices 
                cnt = cnt + 1 
                this%face2vertices(kk,jj) = this%IDNO(cnt) 
            end do 
            if(cnt_vertices < max_vertices) then 
                this%face2vertices(cnt_vertices+1:max_vertices,jj) = -1 
            end if
        end do
            
    end subroutine 

    subroutine make_face2cells(this) 
        type(scf_grid_t) :: this 
        integer :: jj, jjmx 

        jjmx = this%NFACE 

        allocate(this%face2cells(2,jjmx)) 
        !$omp parallel do
        do jj = 1,jjmx 
            this%face2cells(1,jj) = this%IE1(jj) 
            this%face2cells(2,jj) = this%IE2(jj)    !-1のときは存在しない 
        end do 
        !$omp end parallel do
    end subroutine 

    subroutine make_cell2faces(this) 
        type(scf_grid_t) :: this 
        integer :: ii, iimx, jj, jjmx, cnt , cnt_ref
        logical,allocatable :: is_faces(:)

        iimx = this%NELEM 
        jjmx = this%NFACE 

        allocate(is_faces(jjmx)) 
        allocate(this%cell2faces(50,iimx))
        cnt_ref = 0
        do ii = 1,iimx
            ! print*,ii ,'/',iimx
            cnt = 0 
            is_faces(:) = .false. 
            !$omp parallel do 
            do jj = 1,jjmx
                if(this%face2cells(1,jj) == ii-1) then      !要素番号にするために-1
                    is_faces(jj) = .true. 
                else if(this%face2cells(2,jj) == ii-1) then 
                    is_faces(jj) = .true. 
                end if 
            end do 
            !$omp end parallel do

            cnt = count(is_faces) 
            if(cnt > cnt_ref) cnt_ref = cnt 

            cnt = 0 
            do jj = 1, jjmx 
                if(is_faces(jj)) then 
                    cnt = cnt + 1 
                    this%cell2faces(cnt,ii) = jj-1      !面番号にするために-1
                end if 
            end do
            if(cnt < cnt_ref) then 
                this%cell2faces(cnt+1:50,ii) = -1 
            end if 
           
        end do
        ! print*,'cnt_ref =',cnt_ref,'(<50)'

    end subroutine

    subroutine output_txt(this,step)  
        type(scf_grid_t) :: this 
        integer :: step 
        
        !!必要ないものはコメントアウト!!
        if (step == 0) then 
            call output_name_list(this) 
            call output_grid_information(this) 
        end if 
        
        !call output_VOF(this,step) 

    end subroutine


    subroutine output_name_list(this)
        type(scf_grid_t) :: this 
        integer :: i, unit 

        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/namelist.txt',status='replace')
            write(unit,'(A,I8)') '接点数:',this%NODES
            write(unit,'(A,I8)') '要素境界面数:',this%NFACE
            write(unit,'(A,I8)') '要素数:',this%NELEM
            write(unit,'(A)')
            write(unit,'(A)')'EC_Sclar_Data' 
            do i = 1,this%EC_scalar_data_count
                write(unit,'(2A)') this%EC_Scalars(i)%name , this%EC_Scalars(i)%abbreviated_name 
            end do 
            write(unit,'(A)')
            write(unit,'(A)')'EC_Vectors_Data' 
            do i = 1,this%EC_vector_data_count
                write(unit,'(2A)') this%EC_Vectors(i)%name , this%EC_Vectors(i)%abbreviated_name 
            end do 
            write(unit,'(A)') 
            write(unit,'(A)')'FC_Sclar_Data' 
            do i = 1,this%FC_scalar_data_count
                write(unit,'(2A)') this%FC_Scalars(i)%name , this%FC_Scalars(i)%abbreviated_name 
            end do 
            write(unit,'(A)') 
            write(unit,'(A)')'FC_Vectors_Data' 
            do i = 1,this%FC_vector_data_count
                write(unit,'(2A)') this%FC_Vectors(i)%name , this%FC_Vectors(i)%abbreviated_name 
            end do 
        close(unit) 
        print*,'output namelist.txt'
    end subroutine


    subroutine output_grid_information(this) 
        type(scf_grid_t) :: this 
        integer :: i, unit 

        !セル重心の書き出し
        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/cell_centers.txt',status='replace')
            write(unit,'(A,I8)') '要素数:',this%NELEM
            write(unit,'(A)') '要素重心座標: X, Y, Z'
            !要素の中心座標 
            do i = 1, this%NELEM
                write(unit,'(3F12.5)') this%CCE_X(i),this%CCE_Y(i),this%CCE_Z(i) 
            end do 
        close(unit)
        print*,'output cell_centers.txt'
    

        !節点座標の書き出し
        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/nodes.txt',status='replace')
            write(unit,'(A,I8)') '節点数:',this%NODES
            write(unit,'(A)') '節点座標: X, Y, Z'
            !節点の座標 
            do i = 1, this%NODES
                write(unit,'(3F12.5)') this%CAN_X(i),this%CAN_Y(i),this%CAN_Z(i) 
            end do 
        close(unit)
        print*,'output nodes.txt'

        !面情報の書き出し
        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/IE1&IE2.txt',status='replace')
            write(unit,'(A,I8)')'要素界面数:',this%NFACE 
            write(unit,'(A)') 'IE1(境界面裏側の要素番号),IE2(境界面表側の要素番号）(要素番号0〜) ' 
            do i = 1, this%NFACE 
                write(unit,'(3I8)') this%IE1(i),this%IE2(i) 
            end do 
        close(unit)
        print*,'output IE1&IE2.txt'


        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/NDNUM.txt',status='replace')
            write(unit,'(A,I8)')'要素界面数:',this%NFACE 
            write(unit,'(A)') 'NDNUM(界面を構成する節点数)' 
            do i = 1, this%NFACE 
                write(unit,'(I8)') this%NDNUM(i) 
            end do 
        close(unit)
        print*,'output NDNUM.txt'

        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/IDNO.txt',status='replace')
            write(unit,'(A,I8)')'要素界面数:',this%NFACE 
            write(unit,'(A,I9)') 'NDTOT(全界面を構成する節点数):',this%NDTOT  
            write(unit,'(A)') '界面を構成する節点番号（IE1からIE2へ向かう方向へ右ねじ周り）' 
            do i = 1, this%NDTOT
                write(unit,'(I8)') this%IDNO(i)
            end do 
        close(unit) 
        print*,'output IDNO.txt'

        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/face2vertices.txt',status='replace')
            write(unit,'(A,I8)')'要素界面数:',this%NFACE 
            write(unit,'(A)') '界面を構成する節点番号（-1は存在しないことを表す）'
            do i = 1,this%NFACE 
                write(unit,'(10I8)') this%face2vertices(:,i) 
            end do 
        close(unit) 
        print*,'output face2vertices.txt' 
    end subroutine


    subroutine output_VOF(this,step)  
        type(scf_grid_t) :: this 
        integer :: i, unit 
        integer :: step 
        character(len=13) :: fn 

        !VOF値の書き出し
        write(fn,'("VOF_",i5.5,".txt")')step
        open(newunit=unit,file='./'//trim(this%case_name)//'/txt/VOF/'//fn ,status='replace')
            write(unit,'(A,I8)') '要素数:',size(this%EC_Scalars(2)%data)
            write(unit,'(A)') this%EC_Scalars(2)%name 
            !要素のVOF値 
            do i = 1, size(this%EC_Scalars(2)%data)
                write(unit,'(3F12.5)') this%EC_Scalars(2)%data(i)  
            end do 
        close(unit)
        print*,'output',fn  
    end subroutine

end module
