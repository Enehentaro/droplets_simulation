# How to use memo:
# $ cd build
# $ cmake .. -G "MinGW Makefiles" -D CMAKE_Fortran_COMPILER=ifort -D CMAKE_BUILD_TYPE=debug

cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)

project(droplet Fortran)

# #ifort
# if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#     set(CMAKE_Fortran_FLAGS         "-O3 -p -mcmodel=medium -shared-intel")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -fpe0 -g -warn all -check all -standf95 -traceback")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
# endif()
# #gfortran
# if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
#     set(CMAKE_Fortran_FLAGS         "-O3 -mcmodel=medium")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace") #  -p -Wall
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
# endif()

#コンパイルオプション指定用関数
function(SETCOMPILEOPTIONS name)

    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback>
            $<$<CONFIG:Release>:-O2 -fPIC>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        target_compile_options(${name} PRIVATE
            # $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace>
            $<$<CONFIG:Debug>: -O0 -fcheck=all -fbacktrace>
            $<$<CONFIG:Release>: -O2>
        )
    endif()

endfunction(SETCOMPILEOPTIONS)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

enable_testing()

add_subdirectory(lib)
add_subdirectory(src)

add_subdirectory(test)
