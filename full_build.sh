#!/bin/sh

# Based on comments in PluginVsInstall.cmake, 2.10.2
# was compiled with /usr/bin/{gcc,g++} and
# 2.13.0 was compiled with /usr/bin/clang{,++}

# #Linux
# BIN=$HOME/software/visit2_10_3_local/bin/
BIN=$HOME/software/visit3.2.1/bin

rm -r CMake* Makefile cmake_install.cmake

$BIN/xml2info -clobber ugrid.xml
$BIN/xml2cmake -clobber ugrid.xml

# the stdlib flag forces clang++ / g++ to use GNU headers on libraries,
# which have namespaces compatible with the precompiled binaries for VisIt.

# This is related to clang using libc++ (LLVM, with __1:: inline
# namespace) vs. libstdc++ (GNU, no inline namespace)
# cmake -DCMAKE_CXX_FLAGS="-stdlib=libstdc++" . && make VERBOSE=1
cmake . && make VERBOSE=1


