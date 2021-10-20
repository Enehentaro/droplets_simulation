module vector_m
    implicit none
    
    contains

    function cross_product(a, b) result(cross)
        real,intent(in) :: a(3), b(3)
        real cross(3)

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function

    ! function dot_product(a, b) result(inner)
    !     real,intent(in) :: a(3), b(3)
    !     real inner

    !     inner = sum(a(:)*b(:))

    ! end function
end module vector_m