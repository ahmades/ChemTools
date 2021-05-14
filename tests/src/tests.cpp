#define CATCH_CONFIG_RUNNER

#include "config.hpp"

#include <catch2/catch.hpp>

int main(int argc, char* argv[]) {
    
  Catch::Session().run(argc, argv);
  
  return EXIT_SUCCESS;
}
