# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.21)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(K100SIM)

set(CMAKE_BUILD_TYPE None)
#set(WITH_GEANT4_UIVIS OFF)
message("cmake install prefix: ${CMAKE_INSTALL_PREFIX}")
message("source directory: ${CMAKE_CURRENT_SOURCE_DIR}")
set(G4WORKDIR "$ENV{G4WORKDIR}")
message("prefix directory: ${G4WORKDIR}")
execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND git describe --abbrev=100 --always OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
message("git version: ${GIT_VERSION}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__GIT_VERSION=\\\"${GIT_VERSION}\\\"")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__GIT_VERSION=\\\"${GIT_VERSION}\\\"")

#----------------------------------------------------------------------------
#Find HDF5

set (LIB_TYPE STATIC) # or SHARED
string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

get_directory_property(_vars_before VARIABLES)
find_package (HDF5 NAMES hdf5 COMPONENTS C ${SEARCH_TYPE})
get_directory_property(_vars VARIABLES)
# find_package (HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

list(REMOVE_ITEM _vars _vars_before ${_vars_before})
foreach(_var IN LISTS _vars)
    message(STATUS "${_var} = ${${_var}}")
endforeach()
#----------------------------------------------------------------------------
#Setup the HDF5 include directories 

include_directories(${HDF5_INCLUDE_DIR})
set (LINK_LIBS ${LINK_LIBS} ${HDF5_C_${LIB_TYPE}_LIBRARY})

#----------------------------------------------------------------------------
#Find ROOT

#list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}") //only if ROOT made with CMake
list(APPEND CMAKE_MODULE_PATH "$ENV{ROOTSYS}/etc/cmake/")
list(APPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
message(STATUS "CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH})
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_PREFIX_PATH})
get_directory_property(_vars_before VARIABLES)
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})
get_directory_property(_vars VARIABLES)


if (ROOT_FOUND)
    message(STATUS "Libraries for ROOT found")
    message(STATUS "Package Version: " ${ROOT_VERSION})
    message(STATUS "Include DIRS: " ${ROOT_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Failed to find ROOT.")
endif()

#list(REMOVE_ITEM _vars _vars_before ${_vars_before})
#foreach(_var IN LISTS _vars)
#    message(STATUS "${_var} = ${${_var}}")
#endforeach()

#----------------------------------------------------------------------------
#Setup the ROOT include directories 

#include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  message("using visualization")
  get_directory_property(_vars_before VARIABLES)
  find_package(Geant4 REQUIRED ui_all vis_all)
  get_directory_property(_vars VARIABLES)
else()
  message("NOT using visualization")
  get_directory_property(_vars_before VARIABLES)
  find_package(Geant4 REQUIRED)
  get_directory_property(_vars VARIABLES)
endif()

if (Geant4_FOUND)
    message(STATUS "Libraries for Geant4 found")
    message(STATUS "Input directory: " ${Geant4_DIR})
    message(STATUS "Package Version: " ${Geant4_VERSION})
    message(STATUS "Include DIRS: " ${Geant4_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Failed to find Geant4.")
endif()


list(REMOVE_ITEM _vars _vars_before ${_vars_before})
foreach(_var IN LISTS _vars)
    message(STATUS "${_var} = ${${_var}}")
endforeach()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(k100Sim k100Sim.cc ${sources} ${headers})
target_link_libraries(k100Sim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(K100SIM_SCRIPTS
  macros/srcCloseup.mac
  macros/srcRayTrace.mac
  macros/pubeTest.mac
  macros/pube_coin.mac
  sources/pube.dat
  sources/pube_skim.dat
  )

foreach(_script ${K100SIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(K100SIM DEPENDS k100Sim)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS k100Sim 
	RUNTIME DESTINATION ${G4WORKDIR}/bin/Linux-g++)

#----------------------------------------------------------------------------
# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
