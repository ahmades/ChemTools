#include <string>
#include <memory>
#include <map>
#include <stdexcept>
#include <optional>

#include <cassert>

#include "config.hpp"

#include "chemistry.hpp"

namespace Chemistry {
  
  void CleanUp() {
    Cantera::appdelete();
  }

  //
  // -------------- chemisty interface
  //

  IChemistry::IChemistry()
    : thermo( std::make_unique<Cantera::ThermoPhase>() )
    , trans( std::make_unique<Cantera::Transport>() )
    , kinetics( std::make_unique<Cantera::Kinetics>() )
  {}

  Cantera::ThermoPhase* IChemistry::ThermoPtr() {
    return thermo.get();
  }
  
  Cantera::Transport* IChemistry::TransPtr() {
    return trans.get();
  }

  Cantera::Kinetics* IChemistry::KineticsPtr() {
    return kinetics.get();
  }

  std::unique_ptr<Cantera::Transport>& IChemistry::TransUPtr() {
    return trans;
  }

  std::unique_ptr<Cantera::Kinetics>& IChemistry::KineticsUPtr() {
    return kinetics;
  }
    
  //
  // -------------- simplest IChemistry instance only has thermo
  //

  Thermo::Thermo( std::string const& mechanism_ )
    : mechanism( mechanism_ )
  {}
  
  void Thermo::Create() {
    thermo.reset( Cantera::newPhase( mechanism ) );
  }
  
  //
  // -------------- decorator
  //
  Decorator::Decorator( std::unique_ptr<IChemistry> chemistry_ )
    : chemistry( std::move( chemistry_ ) )
  {}
  
  void Decorator::Create() {
    chemistry->Create();
  }

  Cantera::ThermoPhase* Decorator::ThermoPtr() {
    return chemistry->ThermoPtr();
  }
  
  Cantera::Transport* Decorator::TransPtr() {
    return chemistry->TransPtr();
  }

  Cantera::Kinetics* Decorator::KineticsPtr() {
    return chemistry->KineticsPtr();
  }
    
  std::unique_ptr<Cantera::Transport>& Decorator::TransUPtr() {
    return chemistry->TransUPtr();
  }

  std::unique_ptr<Cantera::Kinetics>& Decorator::KineticsUPtr() {
    return chemistry->KineticsUPtr();
  }
     
  //
  // -------------- attaches transport
  //
  Transport::Transport( std::unique_ptr<IChemistry> chemistry
                        , TransportModel model_ /*= TransportModel::MixtureAveraged*/ )
    : Decorator( std::move( chemistry ) )
    , model( model_ )
  {}
  
  void Transport::Create() {
    Decorator::Create();
    std::unique_ptr<Cantera::Transport>& trans_uptr = Decorator::TransUPtr();
    trans_uptr.reset( Cantera::newTransportMgr( TransportModelMap.at( model )
                                                , Decorator::ThermoPtr() ) );
    
  }
  
  //
  // -------------- attaches kinetics
  //
  Kinetics::Kinetics(  std::unique_ptr<IChemistry> chemistry )
    : Decorator( std::move( chemistry ) )
  {}
  
  void Kinetics::Create() {
    Decorator::Create();
    Cantera::ThermoPhase* const thermo_ptr = Decorator::ThermoPtr();
    std::vector<Cantera::ThermoPhase*> const phases{ thermo_ptr };
    std::unique_ptr<Cantera::Kinetics>& kinetics_uptr = Decorator::KineticsUPtr();
    kinetics_uptr.reset( Cantera::newKineticsMgr( thermo_ptr->xml(), phases ) );
  }

  //
  // -------------- convenience factory function
  //
  std::unique_ptr<IChemistry> Create( std::string const& mechanism
                                      , Type type
                                      , TransportModel transport_model
                                      /*= TransportModel::MixtureAveraged*/ ) {
    // attach requested components
    std::unique_ptr<IChemistry> chemistry( nullptr );
    std::unique_ptr<Thermo> thermo( std::make_unique<Thermo>( mechanism ) );
    if( thermo ) {
      switch( type ) {
      case Type::Basic :
        chemistry = std::move( thermo );
        break;
      case Type::Transport:
        chemistry = std::make_unique<Transport>( std::move( thermo ), transport_model );
        break;
      case Type::Kinetics:
        chemistry = std::make_unique<Kinetics>( std::move( thermo ) );
        break;
      case Type::TransportAndKinetics:
        chemistry = std::make_unique<Kinetics>
          ( std::make_unique<Transport>( std::move( thermo ), transport_model ) );
        break;
      default:
        return nullptr;
      }
    } else {
      return nullptr;
    }
    
    // create chemistry
    if( chemistry ) {
      chemistry->Create();
      return chemistry;
    } else {
      return nullptr;
    }
  }

} // namespace Chemistry
