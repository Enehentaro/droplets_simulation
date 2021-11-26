!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine dropletManagement
    use dropletGroup_m
    implicit none
    ! type(dropletGroup) dGroup
    ! integer, save :: cnt = 0

    ! cnt = cnt + 1
    ! if(mod(cnt,100) == 0) then    !100ステップおきに飛沫を10個配置
    !     dGroup = generate_dropletGroup(10)
    !     call mainDroplet%append(dGroup)
    ! end if
    
end subroutine dropletManagement