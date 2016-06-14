#!/bin/sh

#OSX
#BIN=/Applications/VisIt.app/Contents/Resources/bin
#export CC=/usr/local/bin/gcc-4.2
#export CXX=/usr/local/bin/g++-4.2

#Linux
BIN=$HOME/software/visit/bin
export CC=gcc
export CXX=g++

rm -r CMake* Makefile cmake_install.cmake

$BIN/xml2info -clobber ugrid.xml
$BIN/xml2cmake -clobber ugrid.xml

cmake . && make



