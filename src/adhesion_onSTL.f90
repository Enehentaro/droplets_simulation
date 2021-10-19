module adhesion_onSTL_m
    use stl_reader
    implicit none

    type, extends(face_inSTL) :: face_t
        real AB(3), BC(3), CA(3)
    end type face_t

    type(face_t), allocatable :: faceShape(:)

    contains

    subroutine read_faceShape(path)
        implicit none
        type(face_inSTL), allocatable :: solid1(:), solid2(:), solid(:)
        character(*), intent(in) :: path
        character(len=99), allocatable :: stl_fname(:)
        integer i, num_face, num_file
    
        stl_fname = read_stl_fName(path)
        num_file = size(stl_fname)
        print*, 'num_stlFile=', num_file
        solid1 = read_stl(trim(path)//trim(stl_fname(1)))
        solid = solid1
        i = 1
        do while(i < num_file)
            if(allocated(solid)) deallocate(solid)

            i = i + 1
            solid2 = read_stl(trim(path)//trim(stl_fname(i)))

            solid = [solid1, solid2]

            deallocate(solid1, solid2)
            solid1 = solid

        end do

        num_face = size(solid)
        allocate(faceShape(num_face))
        faceShape%face_inSTL = solid
        print*, 'Total Faces on STL:', num_face

    end subroutine read_faceShape

    function read_stl_fName(path) result(stl_fnames)
        character(*), intent(in) :: path
        character(len=99), allocatable :: stl_fnames(:)
        integer n_unit, num_rec, ios, i
        character(len=99) str

        num_rec = 0

        open(newunit=n_unit, file=trim(path)//'stl_list.txt', status='old')
            do        !レコード数を調べるループ
                read(n_unit, *, iostat=ios) str !ファイル終端であればiosに-1が返る
                if((trim(str) == '').or.(ios/=0)) exit    !空白か終端であればループ脱出
                num_rec = num_rec + 1    !カウント
            end do

            allocate(stl_fnames(num_rec))

            rewind(n_unit)

            do i = 1, num_rec
                read(n_unit, *) stl_fnames(i) 
            end do

        close(n_unit)

    end function read_stl_fName

    subroutine set_faceShape
        implicit none
        integer i, num_face
        type(face_inSTL) :: face_

        num_face = size(faceShape)

        do i = 1, num_face
            face_ = faceShape(i)%face_inSTL
            faceShape(i)%AB(:) = face_%nodes(2)%coordinate(:) - face_%nodes(1)%coordinate(:)
            faceShape(i)%BC(:) = face_%nodes(3)%coordinate(:) - face_%nodes(2)%coordinate(:)
            faceShape(i)%CA(:) = face_%nodes(1)%coordinate(:) - face_%nodes(3)%coordinate(:)
        end do

    end subroutine set_faceShape

    logical function adhesion_check_onSTL(X)
        use vector_m
        integer i, num_face
        real, intent(in) :: X(3)
        real r_vec(3), inner, AS(3), BS(3), CS(3), Across(3), Bcross(3), Ccross(3)

        adhesion_check_onSTL = .false.

        num_face = size(faceShape)
        !$OMP parallel do private(r_vec, inner, AS, BS, CS, Across, Bcross, Ccross)
        do i = 1, num_face

            r_vec(:) = X(:) - faceShape(i)%nodes(1)%coordinate(:)  !点Aから飛沫への位置ベクトル
            inner = dot_product(r_vec, faceShape(i)%n_vector(:))    !位置ベクトルと法線ベクトルとの内積
            if (abs(inner) > 1.0d-2) cycle
            AS(:) = r_vec(:) - inner * faceShape(i)%n_vector(:)  !位置ベクトルを面へ投影
            BS(:) = - faceShape(i)%AB(:) + AS(:)
            CS(:) = faceShape(i)%CA(:) + AS(:)
            Across = cross_product(faceShape(i)%AB, AS)
            Bcross = cross_product(faceShape(i)%BC, BS)
            Ccross = cross_product(faceShape(i)%CA, CS)

            !三角形面の内部にあるか判定
            if((dot_product(Across, Bcross) > 0.0).and.(dot_product(Across, Ccross) > 0.0)) then
                adhesion_check_onSTL = .true.
            end if

        end do
        !$OMP end parallel do

    end function adhesion_check_onSTL

end module adhesion_onSTL_m