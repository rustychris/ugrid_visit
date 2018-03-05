#!/bin/sh

# Based on comments in PluginVsInstall.cmake, 2.10.2
# was compiled with /usr/bin/{gcc,g++} and
# 2.13.0 was compiled with /usr/bin/clang{,++}

# Adjust for non-standard installation locations
BIN=/Applications/VisIt.app/Contents/Resources/bin

# High Sierra and VisIt 2.10.2:
# setting nothing will automatically pick up
# a clang which does not cause a warning from
# cmake, but winds up with undefined references to ReadAndProcessDirectory.
# This is discussed in issue #1, and relates to the choice of std c++ library in
# newer Mac OS.

export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

# #Linux
# BIN=$HOME/software/visit/bin
# export CC=gcc
# export CXX=g++

rm -r CMake* Makefile cmake_install.cmake

$BIN/xml2info -clobber ugrid.xml
$BIN/xml2cmake -clobber ugrid.xml

# the stdlib flag forces clang++ / g++ to use GNU headers on libraries,
# which have namespaces compatible with the precompiled binaries for VisIt.

# This is related to clang using libc++ (LLVM, with __1:: inline
# namespace) vs. libstdc++ (GNU, no inline namespace)
cmake -DCMAKE_CXX_FLAGS="-stdlib=libstdc++" . && make VERBOSE=1

