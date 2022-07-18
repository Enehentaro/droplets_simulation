program geometry_test
    use geometry_m
    implicit none
    real, parameter :: tetra(3,4) = reshape([0.,0.,0., 1.,0.,0., 0.,1.,0., 0.,0.,1.], shape(tetra))
    integer, parameter :: imax = 10000
    real rand(3,imax), point(3)
    logical plane_judge, tetra_judge
    integer i

    call random_number(rand)

    do i = 1, imax
        point = rand(:,i)

        plane_judge = (point(3) < plane_equation(point(1), point(2)))
        tetra_judge = insideJudgment_tetra(tetra, point(:))

        if((plane_judge .neqv. tetra_judge)) then
            !テトラの斜面の方程式による内外判定と、テトラの一般的な内外判定の結果が一致しなければエラー
            print*, i, point, plane_judge, tetra_judge
            error stop
        end if

    end do

    contains

    function plane_equation(x,y) result(z)
        real, intent(in) :: x, y
        real z
        z = 1. - x - y
    end function
    
end program geometry_test