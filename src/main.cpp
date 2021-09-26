#include<cstdlib>

#include "chemistry.hpp"

int main( int /*argc*/, char** /*argv*/ ) {

  // when done clean up cantera's xml clutter
  std::atexit( Chemistry::CleanUp );
 
  return EXIT_SUCCESS;
  
}
