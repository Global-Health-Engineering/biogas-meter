# How to make it all run

## How to install REFPROP on M1

1. `conda create -n refprop10 python=3.10`
2. `conda activate refprop10`
3. `pip install six`
4. `brew install gcc`
5. `cd` to your `code` directory
6. `git clone --recursive https://github.com/usnistgov/REFPROP-cmake.git)`
7. `cd REFPROP-cmake`
8. Copy the `FORTRAN` directory into the root of the checked out code
9. Open a console in the root of the cloned repository
10. `mkdir build`
11. `cd build`
12. `where gfortran` - copy the resulting path to the `-DCMAKE_FORTRAN_COMPILER` flag in the next step
12. `cmake .. -DCMAKE_FORTRAN_COMPILER=/path/to/gfortran -DREFPROP_ARM64=ON -DREFPROP_X8664=OFF -DCMAKE_BUILD_TYPE=Release`
13. `cmake --build .`
14. copy the built library `librefprop.dylib` to your REFPROP directory