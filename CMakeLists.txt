
# @file FemSC/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

project (FemSC)

add_subdirectory(headers)
add_subdirectory(headers2018)
add_subdirectory("src vazio")
add_subdirectory(src_2018)
add_subdirectory(tinyxml)


#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/headers2018"
	"${PROJECT_SOURCE_DIR}/headers"
	"${PROJECT_SOURCE_DIR}/tinyxml"
)

## inserting main definitions
#specifying the path to FemSC source code
set(FemSCSOURCE ${FemSC_SRCS} ${FemSC_tinyxml} ${FemSC_HEAD})

set(FemSOURCE2018 ${BASIC_TYPES} ${FemSC_tinyxml} ${HEADERS2018})

add_subdirectory(mainprograms)

## Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
endif()

#define se biblioteca Armadillo deve ser usada
option(USING_ARMADILLO "Biblioteca Armadillo deve ser usada neste projeto" OFF)

if(USING_ARMADILLO)
	find_package( Armadillo REQUIRED NO_MODULE)
	include_directories(${ARMADILLO_INCLUDE_DIRS})
	add_definitions(-DARMADILLO)
endif(USING_ARMADILLO)


#define se biblioteca Eigen3 deve ser usada
option(USING_EIGEN3 "Biblioteca Eigen3 deve ser usada neste projeto" OFF)

if(USING_EIGEN3)
	find_package( Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
	add_definitions(-DEIGEN3)
endif(USING_EIGEN3)


if(USING_ARMADILLO)
	link_libraries(${ARMADILLO_LIBRARIES}) 
endif(USING_ARMADILLO)



option(USING_FAD "Whether the FAD package will be linked with the library" OFF)

if(USING_FAD)
	add_subdirectory("Fad")
endif(USING_FAD)
