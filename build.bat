md build
cd build
cmake .. -G "MinGW Makefile" -D CMAKE_BUILD_TYPE=Debug -D use_OpenMP=ON
make
@REM ctest
cd ..
