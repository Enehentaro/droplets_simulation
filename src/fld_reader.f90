!------------------------------------------------------------------------------
! KIT, AFDET
!------------------------------------------------------------------------------
!
! MODULE:  mod_SctfldReader
!
!> @author
!> Author Name} T.Ikeda, Y.Ida
!
! DESCRIPTION: 
!>  Sc/Tetraのメッシュファイルの解読
!   変数名は基本的にはユーザーズガイドのフォーマットに準拠しています
!
! NOTE:
!>  convert='big_endian'はコンパイラ依存なので，なければコマンド指定
! $ gfortran -fconvert=big-endian read_pre.f90
! $ ifort -convert=big_endian -assume byterecl read_pre.f90
!
! REVISION HISTORY:
!   Ver 0.0_29_May_2021 - Initial Version
!   Ver 0.1_30_May_2021 - 一通り読み込みを完了．各手続きのモジュール組み込み - Ikeda
!   Ver 1.0_30_May_2021 - vtkの変換が成功したので取りあえず暫定版
!
! TODO
!>  構造体の作成
!>  VTK化のための手続きを別に作成
!>  全部privateにしてしまったが，どうするべきか？の検討
!!  各TITLE毎に構造体，読み込みをメンバ手続きにする？
!------------------------------------------------------------------------------
module mod_SctFldReader
    implicit none
    
    private
    !タイトル
    character(32) TITLE
    !サイクル数(=1固定なので実は不要)
    !integer(4) LNX
    !先頭レコード(これもサブルーチン内で定義してしまったので不要)
    !integer(4) IBYTE, IRECN, IRETN
    
    !> 座標番号(=0固定)
    integer(4) LCORD
    
    !> 要素数，1要素あたりの節点数の合計
    !! テトラならNELEM*4?
    integer(4),public ::  NELEM, NTTE

    !> 要素タイプ，要素を構成する節点番号
    !! 読み取りの関係で1次元配列．VTKにするにはこちらで整理する必要あり
    integer(4),allocatable,public :: IETYP(:), NDNO(:)

    !> グループ番号
    integer(4),allocatable,public :: GRP(:)

    !> 物性番号
    !! 多分使えないけど念のため...
    integer(4),allocatable :: MAT(:)

    !> 節点総数, 節点座標
    integer(4),public :: NNODS
    real(8),allocatable,public :: CDN_X(:), CDN_Y(:), CDN_Z(:)

    !> 流速各成分
    real(8),allocatable,public :: VEL_X(:), VEL_Y(:), VEL_Z(:)
    integer(4) NDATA

    !> SctRegion用変数
    ! integer(4) NRGN, NLEN
    ! character,allocatable :: LRGN

    !> 総称手続き
    interface get_data_array
        module procedure get_data_array_int32
        module procedure get_data_array_float64
    end interface get_data_array

    public open_readFLD, close_fld, read_Header_data, read_Main_data

contains

!> 読みこみようにファイルをopen. 装置番号の取得
subroutine open_readFLD(unit, filename)
    implicit none
    character(*),intent(in) :: filename
    integer(4),intent(out) :: unit
    integer(4) ios

    open(newunit = unit, file = filename, form = 'unformatted', &
         access='sequential', convert = 'big_endian', action='read', &
         iostat=ios)
        
        if(ios /= 0) then
            print*,'cannot open file', filename
            return
        else
            write(*,"('Opening fluid file(.fld)...')")
        endif
    
end subroutine open_readFLD

!> ファイルのclose. ファイル操作終了のお知らせを表示
subroutine close_fld(unit)
    implicit none
    integer(4),intent(in) :: unit

    write(*,"('End of Operation: close .fld file')")
    close(unit)
    
end subroutine close_fld

!>序文データの読み取り
subroutine read_Header_data(unit)
    implicit none
    integer(4),intent(in) :: unit
    integer(4)  header_num, NNAMS, n
    character(8)  header_text
    character(32) title_text

    read(unit) header_text          ; print *, header_text
    read(unit) header_num           ; print *, header_num
    ! 序文データの読み取り
    ! おそらくvtk化するためには不要なので
    ! 無限ループでHeaderDataEndを検出した時点でexit
    ! TITLEの画面出力は行う
    do 
        read(unit) title_text ; print *, title_text
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
end subroutine read_Header_data

!/////////////////////////////////////////////////////////////////
!以下本文データ用
!/////////////////////////////////////////////////////////////////

!> 整数型データを読みこみ格納
subroutine get_data_int32(unit,retval)
    implicit none
    integer(4),intent(in) :: unit
    integer(4),intent(out):: retval
    integer(4) ibyte, iretn, irecn 

    read(unit) ibyte, iretn, irecn
    if(iretn*irecn/=1) print*, 'ERROR:get_data_int32'
    read(unit) retval 
    
end subroutine get_data_int32

!> 整数型配列の読み込み
subroutine get_data_array_int32(unit, ret_array, ret_array_size)
    implicit none
    integer(4),intent(in) :: unit
    integer(4),allocatable,intent(out) :: ret_array(:)
    integer(4),intent(in) :: ret_array_size
    integer(4) ibyte, iretn, irecn
    integer(4) irec, L, ios
    integer(4) subrecn

    read(unit) ibyte, iretn, irecn

    subrecn = irecn - 1

    if(.not. allocated(ret_array)) allocate(ret_array(iretn*irecn))

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


end subroutine get_data_array_int32

!> 倍精度実数型配列の読み込み
subroutine get_data_array_float64(unit, ret_array, ret_array_size)
    implicit none
    integer(4),intent(in) :: unit
    real(8),allocatable,intent(out) :: ret_array(:)
    integer(4),intent(in) :: ret_array_size
    integer(4) ibyte, iretn, irecn
    integer(4) irec, L, ios
    integer(4) subrecn

    read(unit) ibyte, iretn, irecn

    subrecn = irecn - 1
    
    if(.not. allocated(ret_array)) allocate(ret_array(iretn*irecn))

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
end subroutine get_data_array_float64

!> データを読み飛ばす処理
subroutine ignore_data(unit)
    implicit none
    integer(4),intent(in) :: unit
    integer(4) ibyte, iretn, irecn
    integer(4) irec

    read(unit) ibyte, iretn, irecn  ;print*, ibyte, iretn, irecn
    do irec = 1, irecn
        read(unit)
    end do

end subroutine ignore_data

!> 本文データの読み取り．とにかく全て読み取るようにしている
subroutine read_Main_data(unit)
    implicit none
    integer(4), intent(in) :: unit
    integer(4) N, NTRY
    character(32) main_data_title

    print*,'---- MAIN DATA START ----'
    !> OVerlapStart_nの読み取り
    read(unit) main_data_title ; print *, main_data_title
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
            call get_data_int32(unit, LCORD)
            read(unit) !0, 0, 0

        case('LS_SurfaceGeometryArray')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !NGFAX
            call ignore_data(unit) !LLEN
            call ignore_data(unit) !LRGNS
            call ignore_data(unit) !NBNNS
            call ignore_data(unit) !IPTYP
            call ignore_data(unit) !IPMAT
            call ignore_data(unit) !NTTSS
            call ignore_data(unit) !NDFA
            read(unit) !0, 0, 0

        case('LS_Elements')
            read(unit) !4, 1, 1
            read(unit) !LNX
            call get_data_int32(unit, NELEM) 
            call get_data_array(unit, IETYP, NELEM) 
            call get_data_int32(unit, NTTE)
            call get_data_array(unit, NDNO, NTTE)  
            read(unit) !0, 0, 0

        case('LS_MatOfElements')
            read(unit) !4, 1, 1
            read(unit) !LNX
            call get_data_int32(unit, NELEM)
            call get_data_array(unit, MAT, NELEM)
            read(unit)

        case('LS_Nodes')
            read(unit) !4, 1, 1
            read(unit) !LNX
            call get_data_int32(unit, NNODS)
            call get_data_array(unit, CDN_X, NNODS)
            call get_data_array(unit, CDN_Y, NNODS)
            call get_data_array(unit, CDN_Z, NNODS)
            read(unit)

        case('LS_VolumeGeometryArray')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !NVORG
            call ignore_data(unit) !LLEN
            call ignore_data(unit) !LRGNS
            call ignore_data(unit) !NELES
            call ignore_data(unit) !IELE
            read(unit)

        case('LS_RegionName&Type')
            call ignore_data(unit) !LNX
            call get_data_int32(unit, NTRY)
            call ignore_data(unit) !NLEN
            do N = 1, NTRY
                call ignore_data(unit) !ITRY
                call ignore_data(unit) !MRGN
            end do
            read(unit)

        case('LS_SFile')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !NLEN
            call ignore_data(unit) !TEXT
            read(unit)

        case('LS_Scalar:PRES')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:TEMP')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:TURK')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)
        
        case('LS_Scalar:TEPS')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:EVIS')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:YPLS')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:HTRC')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Scalar:USTR')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            read(unit)

        case('LS_Vector:VEL')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LVCT
            call ignore_data(unit) !LNAM
            call get_data_int32(unit, NDATA)
            call get_data_array(unit, VEL_X, NDATA)
            call get_data_array(unit, VEL_Y, NDATA)
            call get_data_array(unit, VEL_Z, NDATA)
            read(unit)

        case('LS_Vector:HVEC')
            call ignore_data(unit) !LNX
            call ignore_data(unit) !LVCT
            call ignore_data(unit) !LNAM
            call ignore_data(unit) !NDATA
            call ignore_data(unit) !VAR
            call ignore_data(unit) !VAR
            call ignore_data(unit) !VAR
            read(unit)

        case('OverlapEnd')
            exit

        case default
            print*, 'TITLE_ERROR:', TITLE
            stop

    end select

    end do

    print *, '---- End of file operation ----'

end subroutine read_Main_data

end module mod_SctFldReader

!> テスト用のプログラム
! program main
!     use mod_SctfldReader
!     implicit none
!     integer(4) unit, unit_w
!     !
!     !open(newunit = unit_w, file = 'tut.txt', status='replace')
!     call open_readFLD(unit, 'tutorial_200.fld')
!         call read_Header_data(unit)
!         call read_Main_data(unit)
!     call close_fld(unit)
!     !close(unit_w)

    

!  end program main

!tutorial.preを読ませた結果
!tutorialはユーザーズガイドの操作編のデータで，テトラとプリズムから成る
!NELEM = 358290
!NTTE  = 1689354
!NNODS = 114842