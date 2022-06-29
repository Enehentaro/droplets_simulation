add_library(lib
    caseName.f90
    dropletGenerator.f90
    timeKeeper.f90
    SCTfile_reader.f90
    conditionValue.f90
    dropletMotionSimulation.f90
    path_operator.f90
    unstructured_grid.f90
    adjacency_solver.f90
    filename_mod.f90
    plot3d_operator.f90
    vector.f90
    boxCounter.f90
    flow_field.f90
    simpleFile_reader.f90
    virusDroplet.f90
    dropletEquation.f90
    terminalControler.f90
    vtkMesh_operator.f90
    sort.f90
    kdTree.f90
    kdTree_pointer.f90
    array.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    target_compile_options(lib PRIVATE
        $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback>
        $<$<CONFIG:Release>:-O2 -fPIC>
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    target_compile_options(lib PRIVATE
        $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace>
        $<$<CONFIG:Release>: -O2>
    )
endif()