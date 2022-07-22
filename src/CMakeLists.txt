add_executable(main main.f90)
target_link_libraries(main lib)

add_executable(CUBE2USG CUBE2USG.f90)
target_link_libraries(CUBE2USG lib)

add_executable(boxFlowField boxFlowField.f90)
target_link_libraries(boxFlowField lib)

add_executable(dropletCount dropletCount.f90)
target_link_libraries(dropletCount lib)

add_executable(droplet2CSV droplet2CSV.f90)
target_link_libraries(droplet2CSV lib)

add_executable(initial_translate initial_translate.f90)
target_link_libraries(initial_translate lib)

# if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#     target_compile_options(main PRIVATE
#         $<$<CONFIG:Debug>:-O0 -warn all -check all -traceback>
#         $<$<CONFIG:Release>:-O2 -fPIC>
#     )
# elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
#     target_compile_options(main PRIVATE
#         $<$<CONFIG:Debug>: -O0 -Wall -fcheck=all -fbacktrace>
#         $<$<CONFIG:Release>: -O2>
#     )
# endif()
