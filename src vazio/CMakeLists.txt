# @file FemSC/source/CMakeLists.txt  -- CMake file for the source module

file(GLOB sources *.cpp)

file(GLOB locsources TMatrix.cpp)

set(FemSC_SRCS ${sources} PARENT_SCOPE)

set(BASIC_TYPES ${locsources} PARENT_SCOPE)

##FemSC)SRCS -> nome do arquivo source gerado;