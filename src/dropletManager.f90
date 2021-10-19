module CellCounter_m
    implicit none

    type CellCounter_t
        real volume
        integer counter
    end type CellCounter_t

    type(CellCounter_t), allocatable :: CellCounter(:)

    contains

    subroutine calc_CellVolume
        use unstructuredGrid_mod
        use vector_m
        integer i
        real OA(3), OB(3), OC(3), O(3)

        do i = 1, size(CellCounter)
            select case(CELLs(i)%typeName)
                case('tetra')
                    O(:) = NODEs(CELLs(i)%nodeID(1))%coordinate(:)
                    OA(:) = NODEs(CELLs(i)%nodeID(2))%coordinate(:) - O(:)
                    OB(:) = NODEs(CELLs(i)%nodeID(3))%coordinate(:) - O(:)
                    OC(:) = NODEs(CELLs(i)%nodeID(4))%coordinate(:) - O(:)

                    CellCounter(i)%volume = dot_product(cross_product(OA, OB), OC) / 6.0
            end select
        end do
        
    end subroutine calc_CellVolume

    subroutine output_concentration(FNAME)
        use vtkMesh_operator_m
        character(*), intent(in) :: FNAME
        integer IIMX

        IIMX = size(CellCounter)
        if(.not.allocated(cell_array)) call USG2VTK

        cell_array(0 : IIMX-1)%scalar = CellCounter(1:IIMX)%counter / CellCounter(1:IIMX)%volume

        call output_VTK_mesh(FNAME, data='scalar')
    end subroutine output_concentration

    subroutine USG2VTK
        use unstructuredGrid_mod
        use vtkMesh_operator_m
        integer KKMX, IIMX, IIH, II, KK

        KKMX = size(NODEs)
        if(.not.allocated(node_array)) allocate(node_array(0:KKMX-1))
        do KK = 1, KKMX
            node_array(KK-1)%coordinate(:) = NODEs(KK)%coordinate(:)
        end do
        
        IIMX = size(CELLs)
        if(.not.allocated(cell_array)) allocate(cell_array(0:IIMX-1))
        do II = 1, IIMX
            IIH = size(CELLs(II)%nodeID)
            cell_array(II-1)%nodeID = CELLs(II)%nodeID(1:IIH) - 1
            select case(CELLs(II)%typeName)
                case('tetra')
                    cell_array(II-1)%n_TYPE = 10
                case('prism')
                    cell_array(II-1)%n_TYPE = 13
                case('pyrmd')
                    cell_array(II-1)%n_TYPE = 14
            end select
            cell_array(II-1)%vector(:) = CELLs(II)%flowVelocity(:)
        end do

    end subroutine USG2VTK

end module CellCounter_m

!このサブルーチンは、毎ステップ呼ばれる。自由に編集してよい。
subroutine management_droplet
    use drop_motion_mod
    use CellCounter_m
    implicit none
    character(99) fname
    integer i, refCELL
    
    if(mod(n_time, interval)==0) then
        if(.not.allocated(CellCounter)) then
            allocate(CellCounter(size(CELLs)))
            call calc_CellVolume
        end if
        CellCounter(:)%counter = 0
        do i = 1, size(droplets)
            if(droplets(i)%status < 0) cycle
            refCELL = droplets(i)%refCELL%ID
            CellCounter(refCELL)%counter = CellCounter(refCELL)%counter + 1
        end do
        write(fname,'("'//path%DIR//'concent", i8.8, ".vtk")') n_time
        call output_concentration(trim(fname))
    end if

end subroutine management_droplet