#ifndef RESULTS_OBSERVER_HPP
#define RESULTS_OBSERVER_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <tuple>

#include <cassert>

#include <H5Cpp.h>

namespace Results {
  
  // observer interface
  class Observer {
  public:
    virtual ~Observer() = default;
    virtual void Update( class Subject& ) = 0;
  };

  // Meta of registerd variables
  class Meta {
  public:
  
    enum {
        group    = 0
      , set      = 1
      , name     = 2
      , notation = 3
      , unit     = 4
      , value    = 5
    };

    // scalar result type
    using Scalar = std::tuple< std::string    // group
                               , std::string  // set
                               , std::string  // name
                               , std::string  // notation
                               , std::string  // unit
                               , double*>;    // value

    std::string file_name;
    std::vector<Scalar> scalars;

    Meta()
      : file_name()
      , scalars()
    {}

    void SetFileName( std::string const& file_name_ ) {
      file_name = file_name_;
    }
  
    void Register( std::string const& group
                   , std::string const& set
                   , std::string const& name
                   , std::string const& notation
                   , std::string const& unit
                   , double* const scalar ) {
      scalars.push_back( std::make_tuple( group
                                          , set
                                          , name
                                          , notation
                                          , unit
                                          , scalar ) );
    }

  };

  class ObserverData {
  public:
    Observer* observer;
    size_t update_frequency;
    size_t n_updates;
    
    ObserverData()
      : observer( nullptr )
      , update_frequency{ 1 }
      , n_updates{ 0 }
    {}

    ObserverData( Observer& observer_
                  , size_t update_frequency_ )
      : observer( &observer_ )
      , update_frequency{ update_frequency_ }
      , n_updates{ 0 }
    {}
  };
    
  // subject interface
  class Subject { 
  public:

    Meta results_meta;
    
  
    Subject()
      : results_meta()
      , observers_data()
      , update_frequency()
    {}
  
    virtual ~Subject() = default;
 
    void AttachObserver( Observer& observer
                         , size_t const update_frequency_ ) {
      update_frequency.push_back( update_frequency_ );
      observers_data.emplace_back( ObserverData( observer, update_frequency_ ) );
    }
  
    void DetachObserver( Observer& observer ) {
      update_frequency.pop_back(); 
      observers_data.erase( std::remove_if( observers_data.begin()
                                            , observers_data.end()
                                            , [&]( ObserverData const& data )
                                            { return data.observer == &observer; } )
                            , observers_data.end() );
    }
  
    void NotifyObserver() {
      for( auto& observer_data : observers_data ) {
        if( 0 == ( observer_data.n_updates % observer_data.update_frequency ) ) {
          observer_data.observer->Update( *this );
        }
        observer_data.n_updates++;
      }
    }

    size_t AttachedObservers() {
      return observers_data.size();
    }
  
  private:
  
    std::vector<ObserverData> observers_data;
    std::vector<size_t> update_frequency;
  
  };

} // namespace Results

#endif // RESULTS_OBSERVER_HPP
