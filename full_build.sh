#!/bin/sh

BIN=/Applications/VisIt.app/Contents/Resources/bin

rm -r CMake* Makefile cmake_install.cmake

export CC=/usr/local/bin/gcc-4.2
export CXX=/usr/local/bin/g++-4.2

$BIN/xml2info -clobber ugrid.xml
$BIN/xml2cmake -clobber ugrid.xml

cmake . && make



