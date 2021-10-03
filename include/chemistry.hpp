#ifndef CHEMISTRY_HPP
#define CHEMISTRY_HPP

#include <string>
#include <memory>
#include <map>
#include <stdexcept>
#include <optional>
#include <type_traits>

#include <cassert>

#include "config.hpp"

#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/kinetics.h"

namespace Chemistry {

  // type of transport model when transport is active
  enum class TransportModel {
    None 
    , MixtureAveraged 
    , MultiComponent
  };

  // Map of transport type to string
  static std::map<TransportModel, std::string> const TransportModelMap = {
    { TransportModel::None           , "None"  },
    { TransportModel::MixtureAveraged, "Mix"   },
    { TransportModel::MultiComponent , "Multi" }
  };

  // cleans up xml clutter left out by cantera, should be called on exit
  // lib does not offer a btter solution for this
  inline void CleanUp() {
    Cantera::appdelete();
  }

  // Thermo class, non-copyable
  class Thermo {
  private:
    std::unique_ptr<Cantera::ThermoPhase> m_thermo;
    
  public:
    Thermo()
      : m_thermo( nullptr )
    {}
    
    explicit Thermo( std::string const& mechanism )
      : m_thermo( nullptr )
    {
      Create( mechanism );
    }
    
    Thermo( Thermo const& other ) = delete;
    
    Thermo& operator = ( Thermo const& other ) = delete;
    
    Thermo( Thermo&& other ) = default;
    
    Thermo& operator = ( Thermo&& other ) = default;
    
    ~Thermo() = default;

    void Create( std::string const& mechanism ) {
      if( !m_thermo ) {
        m_thermo.reset( Cantera::newPhase( mechanism ) );
        assert( m_thermo );
      }
    }

    Cantera::ThermoPhase& thermo() const {
      return *m_thermo.get();
    }
    
    Cantera::ThermoPhase& operator()() const {
      return *m_thermo.get();
    }
  };

  // Transport class, non-copyable
  // Creates a Transport object given a ThermoPhase
  class Transport {
  private:
    std::unique_ptr<Cantera::Transport> m_trans;
    
  public:
    Transport()
      : m_trans( nullptr )
    {}
    
    explicit Transport( Cantera::ThermoPhase& thermo
                        , TransportModel transport_model
                        = TransportModel::MixtureAveraged )
      : m_trans( nullptr )
    {
      Create( thermo, transport_model );
    }
    
    Transport( Transport const& other ) noexcept = delete;
    
    Transport& operator = ( Transport const& other ) = delete;
    
    Transport( Transport&& other ) = default;
    
    Transport& operator = ( Transport&& other ) = default;
    
    ~Transport() = default;

    void Create( Cantera::ThermoPhase& thermo
                 , TransportModel transport_model
                 = TransportModel::MixtureAveraged ) {
      if( !m_trans ) {
        m_trans.reset( Cantera::newTransportMgr
                       ( TransportModelMap.at( transport_model ), &thermo ) );
        assert( m_trans );
      }
    }
    
    Cantera::Transport& operator()() const {
      return *m_trans.get();
    }
  };

  // Kinerics class, non-copyable
  // Creates a Kinetics object given a ThermoPhase
  class Kinetics {
  private:
    std::unique_ptr<Cantera::Kinetics> m_kinetics;
    
  public:
    Kinetics()
      : m_kinetics( nullptr )
    {}
    
    explicit Kinetics( Cantera::ThermoPhase& thermo )
      : m_kinetics( nullptr )
    {
      Create( thermo );
    }
    
    Kinetics( Kinetics const& other ) = delete;
    
    Kinetics& operator = ( Kinetics const& other ) = delete;
    
    Kinetics( Kinetics&& other ) = default;
    
    Kinetics& operator = ( Kinetics&& other ) = default;
    
    ~Kinetics() = default;
    
    void Create( Cantera::ThermoPhase& thermo ) {
      if( !m_kinetics ) {
        std::vector<Cantera::ThermoPhase*> const phases{ &thermo };
        m_kinetics.reset( Cantera::newKineticsMgr( thermo.xml(), phases ) );
        assert( m_kinetics );
      }
    }
    
    Cantera::Kinetics& operator()() const {
      return *m_kinetics.get();
    }
  };

  // convenient Thermo + Transport facade
  class ThermoTransport {
  public:
    Thermo thermo;
    Transport transport;

    ThermoTransport() = default;
    
    ThermoTransport( std::string const& mechanism
                     , TransportModel transport_model
                     = TransportModel::MixtureAveraged )
      : thermo( mechanism )
      , transport( thermo(), transport_model )
    {}

    void Create( std::string const& mechanism
                 , TransportModel transport_model
                 = TransportModel::MixtureAveraged ) {
      thermo.Create( mechanism );
      transport.Create( thermo(), transport_model );
    }
    
  };

  // convenien Thermo + Kinetics facade
  class ThermoKinetics {
  public:
    Thermo thermo;
    Kinetics kinetics;

    ThermoKinetics() = default;
    
    ThermoKinetics( std::string const& mechanism )
      : thermo( mechanism )
      , kinetics( thermo() )
    {}

    void Create( std::string const& mechanism ) {
      thermo.Create( mechanism );
      kinetics.Create( thermo() );
    }
  };

  // convenien Thermo +  Transport + Kinetics facade
  class ThermoTransportKinetics {
  public:
    Thermo thermo;
    Transport transport;
    Kinetics kinetics;

    ThermoTransportKinetics() = default;
      
    ThermoTransportKinetics( std::string const& mechanism
                             , TransportModel transport_model
                             = TransportModel::MixtureAveraged )
      : thermo( mechanism )
      , transport( thermo(), transport_model )
      , kinetics( thermo() )
    {}

    void Create( std::string const& mechanism
                 , TransportModel transport_model
                 = TransportModel::MixtureAveraged ) {
      thermo.Create( mechanism );
      transport.Create( thermo(), transport_model );
      kinetics.Create( thermo() );
    }
  };
 
} // namespace Chemistry

#endif //  CHEMISTRY_HPP
