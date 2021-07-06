#define CATCH_CONFIG_RUNNER

#include "config.hpp"

#include <catch2/catch.hpp>

#include "embedded_python/helper.hpp"

int main(int argc, char* argv[]) {

  // Instantiate the python interpreter here to make sure initialisation
  // and finalisation are performed once. Otherwise, seg faults can occur
  // if called from within more than one test.
  CppPy::Instance instance;
   
  Catch::Session().run(argc, argv);
  
  return EXIT_SUCCESS;
}
