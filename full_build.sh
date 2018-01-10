#!/bin/sh

# OSX - used to have homebrew gcc here, but 
# based on comments in PluginVsInstall.cmake, 2.10.2
# was compiled with /usr/bin/{gcc,g++} and
# 2.13.0 was compiled with /usr/bin/clang{,++}

# BIN=/Applications/VisIt2.10.app/Contents/Resources/bin
BIN=/Applications/VisIt2.13.app/Contents/Resources/bin
# High Sierra and VisIt 2.10.2:
# settings nothing will automatically pick up
# a clang which does not cause a warning from
# cmake, but winds up with undefined references to ReadAndProcessDirectory.
# using  /usr/bin/gcc and /usr/bin/g++ has the same outcome.

# Setting CMAKE_C_COMPILER and CMAKE_CXX_COMPILER did nothing
#unset CC
#unset CXX

# Using homebrew gcc, I get other undefined symbols:
#   "vtkDataSetAlgorithm::SetInputData(vtkDataSet*)", referenced from:

# export CC=/usr/local/bin/gcc-4.9
# export CXX=/usr/local/bin/g++-4.9

export CC=/usr/bin/gcc
export CXX=/usr/bin/g++


# #Linux
# BIN=$HOME/software/visit/bin
# export CC=gcc
# export CXX=g++

rm -r CMake* Makefile cmake_install.cmake

$BIN/xml2info -clobber ugrid.xml
$BIN/xml2cmake -clobber ugrid.xml

# -std=c++11 is attempt to fix linking issues with c++ std::string on OSX
# Does not work.
# Seems that this is related to clang using libc++ (LLVM, with __1:: inline
# namespace) vs. libstdc++ (GNU, no inline namespace)
cmake -DCMAKE_CXX_FLAGS="-stdlib=libstdc++" . && make VERBOSE=1


##

# Symbols defined by vtkCommonExecutionModel:
# __ZN19vtkDataSetAlgorithm12SetInputDataEP10vtkDataSet # This is an exact match!
# __ZN19vtkDataSetAlgorithm12SetInputDataEP13vtkDataObject
# __ZN19vtkDataSetAlgorithm12SetInputDataEiP10vtkDataSet
# __ZN19vtkDataSetAlgorithm12SetInputDataEiP13vtkDataObject
