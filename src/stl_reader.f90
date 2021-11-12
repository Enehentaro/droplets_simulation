module stl_reader
    implicit none
    private
    type node
        real coordinate(3)
    end type node
    type, public :: face_inSTL
        type(node) nodes(3)
        real n_vector(3)
    end type

    public read_stl

    contains
    
    function read_stl(fName) result(solid)
        character(*), intent(in) :: fName
        type(face_inSTL), allocatable :: solid(:)
        integer n_unit, ios, num_faces, i
        character str*8, str99*99

        num_faces = 0
        open(newunit=n_unit, file=fName, status='old')
            do        !レコード数を調べるループ
                read(n_unit, *, iostat=ios) str !ファイル終端であればiosに-1が返る
                if(ios/=0) exit    !終端であればループ脱出
                if(trim(str)=='facet') num_faces = num_faces + 1    !面数のカウント
            end do

            rewind(n_unit)

            allocate(solid(num_faces))
            print*, fname, num_faces

            i = 1
            do        !読み取りループ
                read(n_unit, '(A)', iostat=ios) str99 !ファイル終端であればiosに-1が返る
                if(ios/=0) exit    !終端であればループ脱出
                read(str99, *) str
                if(trim(str)=='facet') then
                    read(str99, *) str, str, solid(i)%n_vector(:)
                    read(n_unit, '()')
                    read(n_unit, *) str, solid(i)%nodes(1)%coordinate(:)
                    read(n_unit, *) str, solid(i)%nodes(2)%coordinate(:)
                    read(n_unit, *) str, solid(i)%nodes(3)%coordinate(:)
                    i = i + 1
                end if
            end do
        close(n_unit)

    end function read_stl

end module stl_reader