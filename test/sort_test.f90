!>ヒープソートの機能テスト
program sort_test
    use sort_m
    use array_m
    implicit none
    real, parameter :: test1(10) = [-9.0, -1.0, 0.0, 3.0, 5.0, 7.2, 14.4, 99.9, 122.5, 255.0]
    real, parameter :: test2(11) = [-99.0, -9.0, -1.0, 0.0, 3.0, 5.0, 7.2, 14.4, 99.9, 122.5, 255.0]
    real test3(10000)

    call test_sort(test1, sort_mode='heap_sort')
    call test_sort(test1, sort_mode='merge_sort')

    print '("==============================================================")'

    call test_sort(test2, sort_mode='heap_sort')
    call test_sort(test2, sort_mode='merge_sort')

    print '("==============================================================")'

    block
        integer i
        real rand
        !正解配列を乱数で生成
        test3(1) = 0.
        do i = 2, size(test3)
            call random_number(rand)
            test3(i) = test3(i-1) + rand
        end do
    end block

    call test_sort(test3, sort_mode='heap_sort')
    call test_sort(test3, sort_mode='merge_sort')

    contains

    subroutine test_sort(array_correct, sort_mode)
        !!正解配列（ソート済み配列）を引数に取り、それをシャッフルしたのちにソートを行う
        !!正解配列とソート後の配列を比較し、ソートが機能しているかを検証
        real, intent(in) :: array_correct(:)
        character(*), intent(in) :: sort_mode
        real tmp(size(array_correct))
        type(content_t) c_array(size(array_correct)), c_array_sorted(size(array_correct))
        real array_sorted(size(array_correct))
        integer i

        tmp = FisherYates_shuffle(array_correct)
        ! print *, tmp
        c_array = real2content(tmp)

        select case(sort_mode)
        case('heap_sort')
            c_array_sorted = heap_sort(c_array)
        case('merge_sort')
            c_array_sorted = merge_sort(c_array)
        case default
            error stop
        end select

        array_sorted = c_array_sorted(:)%value

        ! ひとつでも違う要素があればテスト失敗
        if(.not.all(array_sorted == array_correct)) then
            do i = 1, size(array_sorted)
                print '(2(i6, x, f20.16, 4x, "|"))', &
                    c_array(i)%originID, c_array(i)%value, &
                    c_array_sorted(i)%originID, c_array_sorted(i)%value
            end do

            error stop

        end if

    end subroutine
    
end program sort_test