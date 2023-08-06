md build
cd build
cmake .. -G "MinGW Makefiles" -D CMAKE_BUILD_TYPE=Debug -D use_OpenMP=OFF
make
@REM ctest
cd ..
