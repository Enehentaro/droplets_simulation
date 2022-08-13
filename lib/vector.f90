module vector_m
    implicit none
    private
    
    interface cross_product
        !!ベクトルの外積を返す
        module procedure cross_product_dble, cross_product_real
    end interface

    interface normalize_vector
        !!単位ベクトルを返す
        module procedure normalize_vector_dble, normalize_vector_real
    end interface

    public cross_product, normalize_vector
    public norm2_squared

    contains

    function cross_product_dble(a, b) result(cross)
        double precision,intent(in) :: a(3), b(3)
        double precision cross(3)

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function

    function cross_product_real(a, b) result(cross)
        real,intent(in) :: a(3), b(3)
        real cross(3)

        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)

    end function

    function normalize_vector_dble(a) result(normalized)
        double precision,intent(in) :: a(3)
        double precision norm, normalized(3)

        norm = norm2(a)

        normalized(:) = a(:) / norm
        
    end function

    function normalize_vector_real(a) result(normalized)
        real,intent(in) :: a(3)
        real norm, normalized(3)

        norm = norm2(a)

        normalized(:) = a(:) / norm
        
    end function

    real function norm2_squared(a)
        !!ベクトルのL2ノルムの2乗を返す
        !!組み込み関数norm2より約5倍速い(かも)
        real,intent(in) :: a(3)

        norm2_squared = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
        
    end function

end module vector_m