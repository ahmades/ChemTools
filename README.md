# ChemTools

ChemTools is a collection of thermochemical tools for certain chemical and combustion applications.

## Status

This project is currently in the early stages of development. All the necessary infrastructure is in place though. The work on the core applications is in progress.

## Description

Ultimately, the following applications will become available:

- Chemical reactor with one or two properties held constant:
  - pressure,
  - volume,
  - temperature and pressure, or
  - temperature and volume.
- Chemical reactor with user-defined Python scripts providing functions for:
  - the volume and its time derivative as a function of time, or
  - the temperature and pressure as a function of time.
- One-dimensional premixed flames.
- Flamelet library generator using:
  - the β-PDF approach, or
  - the Presumed Mapping Function approach.

## Language, build configuration and testing

ChemTools:

- is written in C++14,
- has CMake support, and
- includes automated unit tests.

## How to build

- To build the code:
  - mkdir build && cd build
  - cmake -DLIB_PATH=/path/to/dependencies /path/to//ChemTools/
- To build the tests (in-source possible for now because some tests use relative paths):
  - cd path/to//ChemTools/
  - mkdir build
  - cmake -DLIB_PATH=/path/to/dependencies -DTEST_VERBOSE=ON ../tests
    Notes:
- Option LIB_PATH adds the specified path to CMAKE_PREFIX_PATH. Please refer to the next section for a complete list of dependencies.
- Option TEST_VERBOSE is used to log some tests results to standard output. By default it is set to OFF.

## Dependencies

ChemTools relies on a number of excellent third-party libraries:

- [Cantera](https://github.com/Cantera/cantera) for the computation of thermodynamic and transport properties, and chemical kinetics,
- [Sundials](https://github.com/LLNL/sundials) for the numerical solution of systems of ordinary differential equations,
- [eigen](https://gitlab.com/libeigen/eigen) for numerical linear algebra computations,
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) for input parsing,
- [units](https://github.com/LLNL/units) for the handling of physical units,
- [hdf5](https://github.com/HDFGroup/hdf5) (C++ API) for the writing of structured simulation results,
- [fmt](https://github.com/fmtlib/fmt) for formatting,
- [Catch2](https://github.com/catchorg/Catch2) for testing,
- [spdlog](https://github.com/gabime/spdlog) for logging (coming soon),
- Python C API, and
- various Boost libraries.

## Platforms

Chemtools is developed and tested under Linux. Build tested with GCC versions 7.3.0 and 11.2.0 and CMake version 3.20.2.

## Documentation

Documentation will be made available as soon some basic applications are published.

## License

ChemTools is distributed under the MIT license.
