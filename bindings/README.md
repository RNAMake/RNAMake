# RNAMake bindings
This subfolder contains the code necessary to make the RNAMake python bindings with pybind11

## Requirements
Python >= 3.0 (not sure on this one, but this is all I have tested in on)
C++17 compatible compiler (have run with both g++-10 and clang++10... no testing for MSVC/MinGW)

## Building
Please follow the below steps

+ Go into `CMakeLists.txt` and use the `set()` function to set `RNAMAKE` to /path/to/RNAMake/. **This is critical for proper building and is a common problem when building the bindings on a new platform**

+ Build the `.so` file independently with `cmake -G Ninja && build`. This takes about 4 minutes on a tenth-gen i7, 32 GB RAM, as a starting point for expectations.

+ Once the build has finished successfully, verify that the `.so` works by entering commandline python and entering `import RNAMake`. If there are no errors, you are good to go! If there are errors, this is usually caused by linking errors, with the most common error being "X is not found". To fix, this, make sure that the referenced symbol actually has a defintion in the code. The name typically looks kind of weird, due to [C++'s name mangling](https://en.wikipedia.org/wiki/Name_mangling). To find out the true name of the symbol, use the `abi::_cxx_demangle()` function found in [`<cxxabi.h>`](https://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html).

+ Run `ninja clean && rm -rf CMakeFiles CMakeCache.txt`. You have to remove these existing files to ensure that the build process in `setup.py` functions properly.

+ Run `pip install --user . `. This took about 5 minutes on a tenth-gen i7, 32 GB RAM.

## Troubleshooting
+ You must have pybind11 cloned. Because it is a submodule, you may have to use the following command: `git submodule update --init --recursive`
+ For pybind11 to work, the name of the file must match what the `PYBIND11` macro is given. For instance, 
`PYBIND(RNAMake,m)` in /RNAMake.cpp **will work**,
but `PYBIND(RNAMake,m)` in /rnamake.cpp or `PYBIND(_RNAMake,m)` in /RNAMake.cpp **will NOT**
+ An `m.def()` directive **cannot** return a smart pointer. If it does, you must derference it in the lamda. 
## TODO
