# min cmake version
cmake_minimum_required(VERSION 3.14.0)

# project description
project(
  chem_tools_tests
  VERSION 1.0.0
  DESCRIPTION "ChemTools automated tests"
  LANGUAGES CXX
  )

# default build type
set(default_build_type "Release")

if(TEST_VERBOSE)
  add_definitions(-DTEST_VERBOSE)
endif()

# executable
add_executable(${PROJECT_NAME})

# sources
target_sources(
  ${PROJECT_NAME}
  PRIVATE
  src/main.cpp
  src/chemistry.cpp
  src/sundials_cvode_iterative_SPGMR.cpp
  src/sundials_cvode_direct_dense.cpp
  src/input_parser.cpp
  src/result_observers.cpp
  src/embedded_python.cpp
  ../src/input_parser.cpp
  ../src/sundials/cvode/interface.cpp
  ../src/results/csv_writer.cpp
  ../src/results/hdf5_writer.cpp
  ../src/results/console_logger.cpp
  ../src/embedded_python_functions.cpp
  )

# includes
target_include_directories(
  ${PROJECT_NAME}
  PRIVATE
  include
  ../include
  )

# c++ standard
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)

# extra compiler options
target_compile_options(${PROJECT_NAME} PRIVATE -s -pedantic -Wall -Wextra -Werror)

# path to dependencies
if(NOT DEFINED LIB_PATH)
  message(FATAL_ERROR
    "Path to dependencies must be specified."
    " All dependencies must be available in the same directory."
    " Use -DLIB_PATH=/path/to/dependenices to set this path."
    )
else()
  list(APPEND CMAKE_PREFIX_PATH ${LIB_PATH})
endif()

# -- Dependency: Catch2
find_package(Catch2 REQUIRED)
if(Catch2_FOUND)
  target_link_libraries(${PROJECT_NAME} Catch2::Catch2)
else()
 message(FATAL_ERROR "Dependency Catch2 not found")
endif()

# -- Dependency: fmt
find_package(fmt REQUIRED)
if(fmt_FOUND)
  target_link_libraries(${PROJECT_NAME} fmt::fmt)
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
  target_link_libraries(${PROJECT_NAME} ${CANTERA_LIBRARY})
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
  string(REPLACE "include/cantera" "include" CANTERA_INCLUDE_DIR ${CANTERA_INCLUDE_DIR})
  target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${CANTERA_INCLUDE_DIR})
else()
 message(FATAL_ERROR "Dependency: Cantera includes not found")
endif()

# -- Dependency: Threads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if(Threads_FOUND)
  target_link_libraries(${PROJECT_NAME} Threads::Threads)
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
    target_link_libraries(${PROJECT_NAME} ${SUNDIALS_LIBRARY})
  else()
    message(FATAL_ERROR "Dependency: Sundials component library ${SUNDIALS_COMPONENT} not found")
  endif()
endforeach()

set(SUNDIALS_INCLUDE_DIR ${CMAKE_PREFIX_PATH}/sundials/include)
if(EXISTS ${SUNDIALS_INCLUDE_DIR})
  target_include_directories(${PROJECT_NAME} PRIVATE  ${SUNDIALS_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Dependency: Sundials includes not found")
endif()

# -- Dependency: yaml-cpp
find_package(yaml-cpp REQUIRED)
if(yaml-cpp_FOUND)
  target_link_libraries(${PROJECT_NAME} yaml-cpp)
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
    ${PROJECT_NAME}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    )
else()
 message(FATAL_ERROR "Dependency Boost not found")
endif()   

# -- Dependency: units
find_package(units REQUIRED)
if(units_FOUND)
  target_link_libraries(${PROJECT_NAME} units-static)
else()
 message(FATAL_ERROR "Dependency units not found")
endif()

# -- Dependency: Eigen3
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
  target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
else()
 message(FATAL_ERROR "Dependency Eigen3 not found")
endif()

#-- Dependency: HDF5
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)
if(HDF5_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${HDF5_CXX_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})
else()
  message(FATAL_ERROR "Dependency HDF5 not found")
endif()

#-- Dependency: Python
# Look preferably for version 3. If not found, version 2.
find_package(Python COMPONENTS Development)
if(Python_FOUND)
  target_include_directories(${PROJECT_NAME} PRIVATE ${Python_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${Python_LIBRARIES})
else()
  message(FATAL_ERROR "Dependency Python not found")
endif()


#find_package(PythonLibs REQUIRED)
#target_include_directories(${PROJECT_NAME} PRIVATE ${PYTHON_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
