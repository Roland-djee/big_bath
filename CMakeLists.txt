# CMakeLists.txt for big_bath
# Dr. Roland Guichard 09/03/2015

cmake_minimum_required(VERSION 2.6)

#### Initialize ################################

# Language
enable_language (Fortran)

# Project name
project (Big_Bath)
set_property(TARGET PROPERTY PROJECT_LABEL ${PROJECT_NAME})

# Build type
set(CMAKE_BUILD_TYPE Release)

#### Packages ##################################

#### Library and executables ###################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin/)
set(main-source_SRC 
src/big_bath.f90 
include/types.f90 
include/constants.f90 
include/read.f90
include/write.f90
library/generate.f90
library/interactions.f90
)

add_executable(big_bath.exe ${main-source_SRC})

set(CMAKE_Fortran_FLAGS "-warn all -check all -traceback -prof-gen=srcpos -align array64byte -fast -parallel -par-report0 ${CMAKE_Fortran_FLAGS}")
