#!/bin/sh

# This is a kludge, temporary fix for some cmake chicanery
# see https://visitbugs.ornl.gov/issues/2611
# and brief mention for this plugin:
# https://github.com/rustychris/ugrid_visit/issues/1

# Adjust this based on Visit version and installation location
# The glob here attempts to grab the right version string which is embedded in the path.
# nonstandard install location or an installation with multiple versions will have to
# manually fix this line
DEP_CMAKE=`echo /Applications/VisIt.app/Contents/Resources/*/darwin-x86_64/include/VisItLibraryDependencies.cmake`

# Seems that the OSX binary distribution for 2.13 has even worse contamination
# of dependencies, and this simple sed expression is not sufficient.

# DEP_CMAKE=/Applications/VisIt2.13.app/Contents/Resources/2.13.0/darwin-x86_64/include/VisItLibraryDependencies.cmake
DEP_CMAKE_O=${DEP_CMAKE}.orig

[ -f $DEP_CMAKE_O ] || cp ${DEP_CMAKE} ${DEP_CMAKE_O}

# libz.dylib is standard on OS X, but erroneously tied to hardcoded paths in the
# cmake dependencies
cat ${DEP_CMAKE_O} | sed -e 's@[-_a-z/A-Z0-9\.]*/libz.dylib@libz.dylib@g' > ${DEP_CMAKE}
