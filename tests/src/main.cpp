#define CATCH_CONFIG_RUNNER

#include "config.hpp"

#include <catch2/catch.hpp>

#include "chemistry.hpp"
#include "embedded_python/helper.hpp"

int main(int argc, char** argv)
{

  // clean up cantera's xml clutter
  std::atexit(Chemistry::CleanUp);

  // Instantiate the python interpreter here to make sure initialisation
  // and finalisation are performed once. Otherwise, seg faults can occur
  // if called from within more than one test.
  CppPy::Instance cpp_py_instance;

  Catch::Session().run(argc, argv);

  return EXIT_SUCCESS;
}
