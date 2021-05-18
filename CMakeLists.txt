cmake_minimum_required(VERSION 3.20.2)

# project description
project(
  PROJ
  VERSION 1.0.0
  DESCRIPTION "ChemTools"
  LANGUAGES CXX
  )

# default build type
set(default_build_type "Release")

# extra compiler options
add_compile_options("-lstdc++ -s -pedantic -Wall -Wextra -Werror")

set(CMAKE_VERBOSE_MAKEFILE ON)

file(GLOB SOURCES "src/*.cpp")

add_executable(exe ${SOURCES})

include_directories(../include)

# c++ standard
set_property(TARGET exe PROPERTY CXX_STANDARD 11)

# path to dependencies
if(NOT DEFINED LIB_PATH)
  message(FATAL_ERROR "Path to dependencies must be specified.\
 All dependencies must be available in the same directory.\
 Use -DLIB_PATH=/pthat/to/dependenices to set the path.")
else()
  set(CMAKE_PREFIX_PATH ${LIB_PATH})
endif()
  
# -- Dependency: fmt
find_package(
  fmt
  REQUIRED
  PATHS ${CMAKE_PREFIX_PATH}/fmt
  )
if(fmt_FOUND)
  message(STATUS "Dependency: fmt found")
  target_link_libraries(exe fmt::fmt)
else()
 message(FATAL_ERROR "Dependency fmt not found")
endif()

#-- Dependency: Cantera
find_library(
  CANTERA_LIBRARY
  NAMES cantera
  HINTS ${CMAKE_PREFIX_PATH}/cantera/lib
  NO_DEFAULT_PATH
  )
if(CANTERA_LIBRARY)
  message(STATUS "Dependency: Cantera library found")
  target_link_libraries(exe ${CANTERA_LIBRARY})
else()
 message(FATAL_ERROR "Dependency: Cantera library not found")
endif()

find_path(
  CANTERA_INCLUDE_DIR
  NAMES thermo.h transport.h kinetics.h core.h
  PATHS ${CMAKE_PREFIX_PATH}/cantera/include NO_DEFAULT_PATH
  PATH_SUFFIXES cantera
  )
if(CANTERA_INCLUDE_DIR)
  message(STATUS "Dependency: Cantera includes found")
  string( REPLACE "include/cantera" "include" CANTERA_INCLUDE_DIR ${CANTERA_INCLUDE_DIR} )
  target_include_directories(exe PRIVATE ${CANTERA_INCLUDE_DIR})
else()
 message(FATAL_ERROR "Dependency: Cantera includes not found")
endif()

# -- Dependency: Threads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if(Threads_FOUND)
  message(STATUS "Dependency: Threads found")
  target_link_libraries(exe Threads::Threads)
else()
  message(FATAL_ERROR "Dependency: Threads not found")
endif()

# -- Dependency: Sundials - CVODE
set(SUNDIALS_REQUIRED_COMPONENTS sundials_cvode sundials_nvecserial)
foreach(SUNDIALS_COMPONENT ${SUNDIALS_REQUIRED_COMPONENTS})
  find_library(
    SUNDIALS_LIBRARY
    NAMES "lib${SUNDIALS_COMPONENT}.so"
    HINTS ${CMAKE_PREFIX_PATH}/sundials/lib
    NO_DEFAULT_PATH
    )
  if(SUNDIALS_LIBRARY)
    message(STATUS "Dependency: Sundials library component ${SUNDIALS_COMPONENT} found")
    target_link_libraries(exe ${SUNDIALS_LIBRARY})
  else()
    message(FATAL_ERROR "Dependency: Sundials component library ${SUNDIALS_COMPONENT} not found")
  endif()
endforeach()

set(SUNDIALS_INCLUDE_DIR ${CMAKE_PREFIX_PATH}/sundials/include)
if(EXISTS ${SUNDIALS_INCLUDE_DIR})
  message(STATUS "Dependency: Sundials includes found")
  target_include_directories(exe PRIVATE  ${SUNDIALS_INCLUDE_DIR})
else()
  message(STATUS "Dependency: Sundials includes found")
endif()
