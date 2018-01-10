#!/bin/sh

# This is a kludge, temporary fix for some cmake chicanery
# see https://visitbugs.ornl.gov/issues/2611
# and brief mention for this pluging:
# https://github.com/rustychris/ugrid_visit/issues/1

# Adjust this based on Visit version and installation location
#DEP_CMAKE=/Applications/VisIt2.10.app/Contents/Resources/2.10.2/darwin-x86_64/include/VisItLibraryDependencies.cmake
DEP_CMAKE=/Applications/VisIt2.13.app/Contents/Resources/2.13.0/darwin-x86_64/include/VisItLibraryDependencies.cmake
DEP_CMAKE_O=${DEP_CMAKE}.orig

[ -f $DEP_CMAKE_O ] || cp ${DEP_CMAKE} ${DEP_CMAKE_O}

# libz.dylib is standard on OS X, but erroneously tied to hardcoded paths in the
# cmake dependencies
cat ${DEP_CMAKE_O} | sed -e 's@[-_a-z/A-Z0-9\.]*/libz.dylib@libz.dylib@g' > ${DEP_CMAKE}
