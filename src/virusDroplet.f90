module virusDroplet_m
    use flow_field
    use equation_mod
    implicit none

    type virusDroplet_t
        double precision :: position(3), velocity(3)=0.d0
        double precision radius, radius_min, initialRadius, deadline
        integer :: status=0, adhesBoundID = 0
        type(reference_cell_t) refCELL

        contains

        procedure evaporation
        procedure motionCalculation
        procedure motionCalculation_onCUBE
        procedure adhesion_onBound
        procedure area_check
        procedure stop_droplet

    end type

    contains

    subroutine evaporation(self) !CALCULATE droplet evaporation
        class(virusDroplet_t) self
        double precision radius_n
      
        if (self%radius <= self%radius_min) return  !半径が最小になったものを除く
    
        radius_n = evaporatin_eq(self%radius)
        
        self%radius = max(radius_n, self%radius_min)
      
    end subroutine

    subroutine motionCalculation(self)
        class(virusDroplet_t) self
        double precision velAir(3)

        velAir(:) = CELLs(self%refCELL%ID)%flowVelocity(:)
    
        call solve_motionEquation(self%position(:), self%velocity(:), velAir(:), self%radius)

        call search_refCELL(real(self%position(:)), self%refCELL%ID)
        
    end subroutine

    subroutine motionCalculation_onCUBE(self)
        class(virusDroplet_t) self
        double precision velAir(3)
        type(reference_cell_t) :: RefC

        RefC = self%refCELL

        velAir(:) = get_velocity_f(RefC%nodeID, RefC%ID)

        call solve_motionEquation(self%position(:), self%velocity(:), velAir(:), self%radius)

        call search_refCELL_onCUBE(real(self%position(:)), self%refCELL)
    
    end subroutine
                    
    subroutine adhesion_onBound(self)
        use vector_m
        class(virusDroplet_t) self
        integer JJ, JB, refCELL
        logical adhesion
        double precision :: r_vector(3), inner

        refCELL = self%refCELL%ID
        adhesion = .false.

        do JJ = 1, size(CELLs(refCELL)%boundFaceID)
            JB = CELLs(refCELL)%boundFaceID(JJ)

            r_vector(:) = self%position(:) - BoundFACEs(JB)%center(:)

            inner = dot_product(r_vector(:), BoundFACEs(JB)%normalVector(:))
            !外向き法線ベクトルと位置ベクトルの内積は、平面からの飛び出し量に相当

            if (inner + self%radius > 0.d0) then
                adhesion = .true.               !(飛び出し量+飛沫半径)がゼロ以上なら付着判定
                self%adhesBoundID = JB       !付着した境界面番号
            end if
        end do

        if (adhesion) call self%stop_droplet()

    end subroutine

    subroutine area_check(self)
        class(virusDroplet_t) self
        logical check
        integer L

        check = .false.
        do L = 1, 3
    
            if(self%position(L) < MIN_CDN(L)) then
                self%position(L) = MIN_CDN(L)
                check = .true.
            else if(self%position(L) > MAX_CDN(L)) then
                self%position(L) = MAX_CDN(L)
                check = .true.
            end if

        end do

        if (check) call self%stop_droplet()

    end subroutine

    subroutine stop_droplet(self, status)
        class(virusDroplet_t) self
        integer, optional :: status

        self%velocity(:) = 0.0d0
        if(present(status)) then
            self%status = status
        else
            self%status = 1
        end if

    end subroutine

end module virusDroplet_m