# PLAN
PLanetesimal ANalyzer is designed to find & locate planetesimals and calculate their characteristics in the output data of numerical simulations by code ATHENA.

## Compile & Run
CMake is needed to generate a Makefile and compile this program. Boost library is also required. 

```bash
➜  PLAN $ ls
CMakeLists.txt README.md      src
➜  PLAN $ mkdir build; cd ./build
➜  build $ cmake ..
-- The C compiler identification is AppleClang 8.1.0.8020042
-- The CXX compiler identification is AppleClang 8.1.0.8020042
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
============================================================
Generating Makefile for serial program...
For parallel build, run cmake with -DPARALLEL=ON
============================================================
-- Boost version: 1.59.0
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/rixin/PLAN
➜  build $ make -j 4
Scanning dependencies of target plan
[ 80%] Building CXX object CMakeFiles/plan.dir/src/tree.cpp.o
[ 80%] Building CXX object CMakeFiles/plan.dir/src/global.cpp.o
[ 80%] Building CXX object CMakeFiles/plan.dir/src/analyses.cpp.o
[ 80%] Building CXX object CMakeFiles/plan.dir/src/main.cpp.o
[100%] Linking CXX executable plan
[100%] Built target plan
➜  build $ ls
CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile  plan
```

Specify the option, "PARALLEL", to build parallel program (MPI library is required).

```bash
➜  PLAN $ ls
CMakeLists.txt README.md      src
➜  PLAN $ mkdir build; cd ./build
➜  build $ cmake -DPARALLEL=ON ..
-- The C compiler identification is AppleClang 8.1.0.8020042
-- The CXX compiler identification is AppleClang 8.1.0.8020042
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
============================================================
Generating Makefile for parallel program
============================================================
-- Found MPI_C: /opt/local/lib/openmpi-mp/libmpi.dylib
-- Found MPI_CXX: /opt/local/lib/openmpi-mp/libmpi_cxx.dylib;/opt/local/lib/openmpi-mp/libmpi.dylib
-- Boost version: 1.59.0
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/rixin/PLAN
➜  build $ make -j 4
Scanning dependencies of target plan
[ 40%] Building CXX object CMakeFiles/plan.dir/src/global.cpp.o
[ 40%] Building CXX object CMakeFiles/plan.dir/src/tree.cpp.o
[ 60%] Building CXX object CMakeFiles/plan.dir/src/analyses.cpp.o
[ 80%] Building CXX object CMakeFiles/plan.dir/src/main.cpp.o
[100%] Linking CXX executable plan
[100%] Built target plan
➜  build $ ls
CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile  plan
```

To clean CMake results and cache, just delete the `build` directory.


