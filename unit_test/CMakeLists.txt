# CMakeLists.txt for big_bath unit tests.
# Dr. Roland Guichard 09/03/2015

cmake_minimum_required(VERSION 2.6)

#### Initialize ################################

# Language
enable_language (Fortran)

# Project name
project(unit_test)
set_property(TARGET PROPERTY PROJECT_LABEL ${PROJECT_NAME})

# Build type
set(CMAKE_BUILD_TYPE Release)

#### Packages ##################################

#### Library and executables ###################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)

set(main-source_SRC
../../include/types.f90 
../../include/constants.f90 
../../include/read.f90
../../include/write.f90
../../library/generate.f90
../../library/interactions.f90
src/unit_test.f90
test/test_read.f90
test/test_generate.f90
test/test_interactions.f90)

add_executable(unit_test.exe ${main-source_SRC})
set(CMAKE_Fortran_FLAGS "-prof-gen=srcpos -check all -warn all -traceback ${CMAKE_Fortran_FLAGS}")





