program array2VTK
    use array_m
    use VTK_operator_m
    implicit none
    real, allocatable :: UVW(:,:)
    type(UnstructuredGrid_inVTK) usg
    character(255) array_fname, output_fname
    character(:), allocatable :: path, casepath
    integer i

    path = '/run/user/1000/gvfs/smb-share:server=nas2.local,share=enehen-master2/konishi/output/array/office11_data/'
    casepath = path // '960_246_aout/output/'
    call usg%read(path//'office_usg.vtk')

    do i = 3000, 5125, 25
        write(array_fname, '(A, "field_" i10.10, ".array")') casepath, i
        call read_2dArray_asBinary(array_fname, UVW)
        write(output_fname, '(A, "field_" i10.10, ".vtk")') casepath, i
        call USG%output(output_fname, cellVector=UVW, vectorName='Velocity')
    end do
    
end program array2VTK
