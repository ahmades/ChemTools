#ifndef RESULTS_CONSOLE_LOGGER_HPP
#define RESULTS_CONSOLE_LOGGER_HPP

#include "results_observer.hpp"

namespace Results {

  class ConsoleLogger: public Observer { 
  public:
    
    ConsoleLogger( Subject& subject_
                   , size_t const update_frequency_ = 1 );
  
    ~ConsoleLogger();
    
    // updates values of registered variables
    void Update( Subject& a_subject ) override;
    
  private:
    
    Subject& subject;
  };
  
} // namespace Results


#endif // RESULTS_CONSOLE_LOGGER_HPP
