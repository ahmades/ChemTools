# min cmake version
cmake_minimum_required(VERSION 3.14.0)

# project description
project(
  PROJ
  VERSION 1.0.0
  DESCRIPTION "ChemTools"
  LANGUAGES CXX
  )

# default build type
set(default_build_type "Release")

# executable
add_executable(exe)

# sources
target_sources(
  exe
  PRIVATE
  src/main.cpp
  src/input_parser.cpp
  src/sundials_cvode.cpp
  )

# includes
target_include_directories(
  exe
  PRIVATE
  include
  )

# c++ standard
set_property(TARGET exe PROPERTY CXX_STANDARD 14)

# extra compiler options
target_compile_options(test_exe PRIVATE -pedantic -Wall -Wextra -Werror)

# path to dependencies
if(NOT DEFINED LIB_PATH)
  message(FATAL_ERROR
    "Path to dependencies must be specified."
    " All dependencies must be available in the same directory."
    " Use -DLIB_PATH=/path/to/dependenices to set this path."
    )
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
  string( REPLACE "include/cantera" "include" CANTERA_INCLUDE_DIR ${CANTERA_INCLUDE_DIR} )
  target_include_directories(exe SYSTEM PRIVATE ${CANTERA_INCLUDE_DIR})
else()
 message(FATAL_ERROR "Dependency: Cantera includes not found")
endif()

# -- Dependency: Threads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if(Threads_FOUND)
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
    target_link_libraries(exe ${SUNDIALS_LIBRARY})
  else()
    message(FATAL_ERROR "Dependency: Sundials component library ${SUNDIALS_COMPONENT} not found")
  endif()
endforeach()

set(SUNDIALS_INCLUDE_DIR ${CMAKE_PREFIX_PATH}/sundials/include)
if(EXISTS ${SUNDIALS_INCLUDE_DIR})
  target_include_directories(exe PRIVATE  ${SUNDIALS_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Dependency: Sundials includes not found")
endif()

# -- Dependency: yaml-cpp
find_package(
  yaml-cpp REQUIRED
  PATHS ${CMAKE_PREFIX_PATH}/yaml-cpp
  )
if(yaml-cpp_FOUND)
  target_link_libraries(exe yaml-cpp)
else()
 message(FATAL_ERROR "Dependency yaml-cpp not found")
endif()

# -- Dependency: Boost components
find_package(
  Boost
  COMPONENTS system filesystem
  REQUIRED
  )
if(Boost_FOUND)
  target_link_libraries(
    exe
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    )
else()
 message(FATAL_ERROR "Dependency Boost not found")
endif()   

# -- Dependency: units
find_package(
  units REQUIRED
  PATHS ${CMAKE_PREFIX_PATH}/units
  )
if(units_FOUND)
  target_link_libraries(exe units-static)
else()
 message(FATAL_ERROR "Dependency units not found")
endif()

# -- Dependency: Eigen3
find_package(
  Eigen3 REQUIRED
  PATHS ${CMAKE_PREFIX_PATH}/eigen
  )
if(Eigen3_FOUND)
  target_link_libraries(test_exe Eigen3::Eigen)
else()
 message(FATAL_ERROR "Dependency Eigen3 not found")
endif()
