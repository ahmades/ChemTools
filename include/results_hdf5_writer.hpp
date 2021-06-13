#ifndef RESULTS_HDF5_WRITER_HPP
#define RESULTS_HDF5_WRITER_HPP

#include <unordered_map>

#include <H5Cpp.h>

#include "results_observer.hpp"

namespace Results {

  class HDF5Writer: public Observer {
  public:
    
    HDF5Writer( Subject& subject_
                , std::string const& path
                , size_t const update_frequency_ = 1 );
    
    ~HDF5Writer();

    // Modifies the default compression level
    void SetCompressionLevel( int const compression_level_ );
    
    // updates the values of registered variables
    void Update( Subject& a_subject ) override;
    
  private:
    
    Subject& subject;
    std::string path;
    std::unique_ptr<H5::H5File> file;
    std::unordered_map< std::unique_ptr<H5::DataSet>, double* > map_scalars;
    int compression_level;

    // bundles the different initialisation steps
    void Initialise();

    // opens the results file
    void OpenFile();

    // creates gropus of registered variables
    void CreateGroups();

    // creates scalar datasets in groups
    void CreateScalarDataSets();
    
  };
  
}

#endif // RESULTS_HDF5_WRITER_HPP
