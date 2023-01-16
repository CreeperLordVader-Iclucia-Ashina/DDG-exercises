# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/deps/geometry-central/deps/eigen-src"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/deps/geometry-central/deps/eigen-build"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/tmp"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/src"
  "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/creeper/ddg-exercises/projects/simplicial-complex-operators/cmake-build-debug/_deps/eigen-subbuild/eigen-populate-prefix/src/eigen-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
