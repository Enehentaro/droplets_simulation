module array_m
    implicit none
    private

    interface isEqual
        module procedure isEqual_int, isEqual_real
    end interface

    public isEqual
    public FisherYates_shuffle

    contains

    !配列が等しいかどうかを判定する関数
    function isEqual_int(a, b) result(isEqual_)
        integer, intent(in) :: a(:), b(:)
        logical isEqual_
        integer i

        isEqual_ = .true.

        do i = 1, size(a)
            if(a(i)/=b(i)) then
                isEqual_ = .false.
                return
            end if
        end do

    end function

    !配列が等しいかどうかを判定する関数
    function isEqual_real(a, b) result(isEqual_)
        real, intent(in) :: a(:), b(:)
        logical isEqual_
        integer i

        isEqual_ = .true.

        do i = 1, size(a)
            if(a(i)/=b(i)) then
                isEqual_ = .false.
                return
            end if
        end do

    end function

    function FisherYates_shuffle(a) result(b)
        real, intent(in) :: a(:)
        real b(size(a)), rand, tmp
        integer i, index

        b = a

        do i = size(b), 2, -1
            call random_number(rand)
            index = int(rand * (i-1)) + 1
            ! print *, index, i

            !SWAP
            tmp = b(index)
            b(index) = b(i)
            b(i) = tmp

        end do

    end function

end module array_m