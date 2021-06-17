#ifndef RESULTS_CSV_WRITER_HPP
#define RESULTS_CSV_WRITER_HPP

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <memory>
#include <iomanip>

#include <cassert>

#include "results/observer.hpp"

namespace Results {

  class CSVWriter: public Observer {
  public:
    
    CSVWriter( Subject& subject_
               , std::string const& path
               , size_t const update_frequency_ = 1 );
    
    ~CSVWriter();
    
    // updates values of registered variables as of line 4
    void Update( Subject& a_subject ) override;
  
  private:
  
    Subject& subject;
    std::string path;
    std::ofstream stream;

    // opens stream
    void OpenFile();
    
    /* inserts the header:  
        - line 1: group.set.names 
        - line 2 : notation
        - line 3: units
     */
    void InsertHeader();
  };
  
} // namespace Results

#endif // RESULTS_CSV_WRITER_HPP
