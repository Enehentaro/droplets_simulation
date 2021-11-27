!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use dropletGroup_m
    implicit none
    ! type(dropletGroup) dGroup

    ! if(mod(n_time, 100) == 0) then    !100ステップおきに飛沫を10個発生させるサンプルコード
    !     dGroup = generate_dropletGroup(10)
    !     call mainDroplet%append(dGroup)
    ! end if
    
end subroutine dropletManagement