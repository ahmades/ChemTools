#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>

#include <fmt/format.h>

#include "results/csv_writer.hpp"

namespace Results {

  CSVWriter::CSVWriter( Subject& subject_
                        , std::string const& path_ )
    : subject( subject_ )
    , path( path_ )
    , stream()
  {
    // attach
    subject.AttachObserver( *this );
    
    // open the file
    OpenFile();
    
    //insert header of registered data
    InsertHeader();
  }
  
  CSVWriter::~CSVWriter() {
    subject.DetachObserver( *this );
  }
  
  void CSVWriter::Update( Subject& a_subject ) {
    if( &a_subject == &subject ) {
      std::vector<Meta::Scalar> const& scalars = subject.results_meta.scalars;
      for( auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars ) {
        stream << *std::get<Meta::value>( *it_scalars );
        if( it_scalars != scalars.end()-1 ) stream << ',';
      }
      stream << '\n';
    }
  }

  void CSVWriter::OpenFile() {
    // check path
    boost::filesystem::path const fs_path( path );
    boost::system::error_code error_code;
    if( !boost::filesystem::exists( fs_path, error_code ) ) {
      BOOST_FILESYSTEM_THROW( boost::filesystem::filesystem_error( "CSVWriter error"
                                                                   , fs_path
                                                                   , error_code ) );
    }
    
    // still need to get absolute path if relative is supplied
    
    // open stream for writing
    stream.open( path + '/' + subject.results_meta.file_name + ".csv"
                 , std::fstream::out
                 | std::fstream::trunc );
  }
  
  void CSVWriter::InsertHeader() {
    std::vector<Meta::Scalar> const& scalars = subject.results_meta.scalars;   
    // names
    for( auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars ) {
      stream << fmt::format( "{}.{}.{}"
                             , std::get<Meta::group>( *it_scalars )
                             , std::get<Meta::set>( *it_scalars )
                             , std::get<Meta::name>( *it_scalars ) );
      
      if( it_scalars != scalars.end()-1 ) stream << ',';
    }
    stream << '\n';

    // notation
    for( auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars ) {
      stream << std::get<Meta::notation>( *it_scalars );
      if (it_scalars != scalars.end()-1 ) stream << ',';
    }
    stream << '\n';
    
    // units
    for( auto it_scalars = scalars.begin(); it_scalars != scalars.end(); ++it_scalars ) {
      stream << std::get<Meta::unit>( *it_scalars );
      if (it_scalars != scalars.end()-1 ) stream << ',';
    }
    stream << '\n';
  }

} // namspace Results
