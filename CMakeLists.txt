cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)

project(droplet Fortran)

#ifort
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-O3 -p -mcmodel=medium -shared-intel")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -fpe0 -g -warn all -check all -standf95 -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
endif()
#gfortran
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "-O3 -p -mcmodel=medium")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
endif()

# $cmake -D CMAKE_Fortran_COMPILER=ifort -D CMAKE_BUILD_TYPE=debug 

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/fortran-modules)

# add_subdirectory(util)
add_subdirectory(src)
# add_subdirectory(lib)
