module geometry_m
    use vector_m
    implicit none
    private

    real, parameter :: delta = 1.e-6
    real, parameter :: frac_1_6 = 1./6.

    public volume_tetra, insideJudgment_tetra, insideJudgment_tetra_check
    
contains

    !>テトラの体積計算。
    function volume_tetra(vertices) result(volume)
        real, intent(in) :: vertices(3,4)
        real volume
        real, dimension(3) :: a,b,c

        a = vertices(:,2) - vertices(:,1)
        b = vertices(:,3) - vertices(:,1)
        c = vertices(:,4) - vertices(:,1)

        volume = frac_1_6 * abs(dot_product(cross_product(a, b), c))

    end function

    !>任意の点がテトラの内部にあるかどうかを判定する。
    !>点でテトラを分割したそれぞれの体積の和が、元々のテトラの体積を上回れば、点はテトラ外部にある。
    !>https://matcha-choco010.net/2018/03/14/point-in-tetrahedron/
    function insideJudgment_tetra(vertices, point) result(isInside)
        real, intent(in) :: vertices(3,4), point(3)
        real volume, vol_sum
        logical isInside

        volume = volume_tetra(vertices)
        vol_sum = partitionedVolume_sum_tetra(vertices, point)

        !点がテトラの内部にあるとき、分割体積和と元々の体積は厳密に一致するはずだが、
        !数値誤差のために少し条件を緩和する
        isInside = (vol_sum <= volume*(1. + delta))

    end function

    !>点でテトラを分割したそれぞれの体積の和を求める
    function partitionedVolume_sum_tetra(vertices, point) result(vol_sum)
        real, intent(in) :: vertices(3,4), point(3)
        real vol_sum, mini_vertices(3,4)
        integer i

        vol_sum = 0.

        do i = 1, 4
            mini_vertices = vertices
            mini_vertices(:,i) = point
            vol_sum = vol_sum + volume_tetra(mini_vertices)
        end do

    end function

    subroutine insideJudgment_tetra_check(vertices, point, vol_sum, volume)
        real, intent(in) :: vertices(3,4), point(3)
        real volume, vol_sum

        volume = volume_tetra(vertices)
        vol_sum = partitionedVolume_sum_tetra(vertices, point)

    end subroutine

end module geometry_m