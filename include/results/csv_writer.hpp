#ifndef RESULTS_CSV_WRITER_HPP
#define RESULTS_CSV_WRITER_HPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <cassert>

#include "results/observer.hpp"

namespace Results
{

  class CSVWriter : public Observer
  {
  public:
    CSVWriter(Subject& subject_, std::string const& path);

    ~CSVWriter();

    // updates values of registered variables as of line 4
    void Update(Subject& a_subject) override;

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
