#ifndef RESULTS_CONSOLE_LOGGER_HPP
#define RESULTS_CONSOLE_LOGGER_HPP

#include "results/observer.hpp"

namespace Results {

  class ConsoleLogger: public Observer { 
  public:
    
    ConsoleLogger( Subject& subject_ );
  
    ~ConsoleLogger();
    
    // updates values of registered variables
    void Update( Subject& a_subject ) override;
    
  private:
    
    Subject& subject;
  };
  
} // namespace Results


#endif // RESULTS_CONSOLE_LOGGER_HPP
