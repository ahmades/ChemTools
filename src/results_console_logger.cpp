#include <fmt/format.h>
#include <fmt/printf.h>

#include "results_console_logger.hpp"

namespace Results {
  
  ConsoleLogger::ConsoleLogger( Subject& subject_
                                , size_t const update_frequency_ )
    : subject( subject_ )
  {
    // attach
    subject.AttachObserver( *this, update_frequency_ );
  }
  
  ConsoleLogger::~ConsoleLogger() {
    subject.DetachObserver( *this );
  }
  
  // logs the values of registered variables
  void ConsoleLogger::Update( Subject& a_subject ) {
    if( &a_subject == &subject ) {
      std::vector<Meta::Scalar> const& scalars = subject.results_meta.scalars;
      
      for( auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars ) {
        fmt::print( "{} =  {:e} {}"
                    , std::get<Meta::notation>( *it_scalars )
                    , *std::get<Meta::value>( *it_scalars )
                    , std::get<Meta::unit>( *it_scalars ) );
        
        if( it_scalars != scalars.end()-1 ) fmt::print( " " );
      }
      fmt::print( "\n" );
    }
  }

} // namsepace Results
