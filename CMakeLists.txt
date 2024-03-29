# How to use memo:
# $ cd build
# $ cmake .. -G "MinGW Makefiles" -D CMAKE_Fortran_COMPILER=ifort -D CMAKE_BUILD_TYPE=debug

cmake_minimum_required(VERSION 3.5)

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

option(use_OpenMP "enable OpenMP directives" OFF)

# パッケージの探索(コンパイラにOpenMPが付属していなければインストール)
if(use_OpenMP)
    find_package(OpenMP REQUIRED)
endif()

#コンパイルオプション指定用関数
function(SETCOMPILEOPTIONS name)

    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>:-O0 -traceback> # -warn all -check all 
            $<$<CONFIG:Release>:-O2 -fPIC>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>: -O0 -g -fcheck=all> # -Wall
            $<$<CONFIG:Release>: -O2>
        )
    endif()

    if(OpenMP_Fortran_FOUND)
        target_compile_options(${name} PRIVATE 
            "${OpenMP_Fortran_FLAGS}"
        )
        target_link_libraries(${name}
            OpenMP::OpenMP_Fortran 
        )
    endif()

endfunction(SETCOMPILEOPTIONS)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

enable_testing()

add_subdirectory(lib)
add_subdirectory(src)

add_subdirectory(test)
