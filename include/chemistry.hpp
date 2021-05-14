#ifndef CHEMISTRY_HPP
#define CHEMISTRY_HPP

#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <stdexcept>
#include <optional>

#include <cassert>

#include "config.hpp"

#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/kinetics.h"

namespace Chemistry {
  
  // type of chemistry to create
  enum class Type {
    Invalid
    , Basic
    , Transport
    , Kinetics
    , TransportAndKinetics
   };

  // type of transport model when transport is active
  enum class TransportModel {
    None 
    , MixtureAveraged 
    , MultiComponent
  };

  // Map of transport type to string
  static std::map<TransportModel, std::string> TransportModelMap = {
    { TransportModel::None           , "None"  },
    { TransportModel::MixtureAveraged, "Mix"   },
    { TransportModel::MultiComponent , "Multi" }
  };

  // cleans up xml clutter left out by cantera, should be called on exit
  // lib does not offer a btter solution for this
  static void CleanUp() {
    Cantera::appdelete();
  }

  //
  // -------------- chemisty interface
  //
  class IChemistry {
  protected:
    std::unique_ptr<Cantera::ThermoPhase> thermo;
    std::unique_ptr<Cantera::Transport> trans;
    std::unique_ptr<Cantera::Kinetics> kinetics;
  
  public:

    IChemistry()
      : thermo( nullptr )
      , trans( nullptr )
      , kinetics( nullptr )
    {};
  
    virtual ~IChemistry() = default;

    virtual void Create() = 0;
 
    virtual std::unique_ptr<Cantera::ThermoPhase>& GetThermo() {
      return thermo;
    }
  
    virtual std::unique_ptr<Cantera::Transport>& GetTrans() {
      return trans;
    }

    virtual std::unique_ptr<Cantera::Kinetics>& GetKinetics() {
      return kinetics;
    }
  };

  //
  // -------------- simplest IChemistry instance only has thermo
  //
  class Thermo: public IChemistry {
  public:

    Thermo( std::string const& mechanism_ )
      : mechanism( mechanism_ )
    {}
  
    ~Thermo() = default;
  
    void Create() override {
      thermo = std::unique_ptr<Cantera::ThermoPhase>( Cantera::newPhase( mechanism ) );
    }
  
  private:
    
    std::string mechanism;
  };

  //
  // -------------- decorator
  //
  class Decorator: public IChemistry {
  public:
  
    Decorator( IChemistry* chemistry_ )
      : chemistry( chemistry_ )
    {}

    ~Decorator() = default;
 
    virtual void Create() {
      chemistry->Create();
    }

    std::unique_ptr<Cantera::ThermoPhase>& GetThermo() override {
      return chemistry->GetThermo();
    }
  
    std::unique_ptr<Cantera::Transport>& GetTrans() override {
      return chemistry->GetTrans();
    }

    std::unique_ptr<Cantera::Kinetics>& GetKinetics() override {
      return chemistry->GetKinetics();
    }
   
  private:
    
    std::unique_ptr<IChemistry> chemistry;

  };
  
  //
  // -------------- attaches transport
  //
  class Transport : public Decorator {
  public:

    Transport( IChemistry* chemistry
               , TransportModel model_ = TransportModel::MixtureAveraged )
      : Decorator( chemistry )
      , model( model_ )
    {}

    ~Transport() = default;
  
    void Create() override {
      Decorator::Create();
      std::unique_ptr<Cantera::Transport>& trans_uptr = Decorator::GetTrans();
      trans_uptr = std::unique_ptr<Cantera::Transport>
        ( Cantera::newTransportMgr( TransportModelMap.at( model )
                                    , Decorator::GetThermo().get() ) );
    }
  
  private:
    
    TransportModel model;
  };

  //
  // -------------- attaches kinetics
  //
  class Kinetics : public Decorator {
  public:
  
    Kinetics( IChemistry* chemistry )
      : Decorator( chemistry )
    {}

    Kinetics() = default;

    void Create() override {
      Decorator::Create();
      std::unique_ptr<Cantera::ThermoPhase> const& thermo_uptr = Decorator::GetThermo();
      std::vector<Cantera::ThermoPhase*> phases{ thermo_uptr.get() };
      std::unique_ptr<Cantera::Kinetics>& kinetics_uptr = Decorator::GetKinetics();
      kinetics_uptr = std::unique_ptr<Cantera::Kinetics>
          ( Cantera::newKineticsMgr( thermo_uptr->xml(), phases ) );
    }
  };

  //
  // -------------- convenience factory
  //
  static std::unique_ptr<IChemistry>
  Create( std::string const& mechanism
          , Type const type
          , TransportModel transport_model = TransportModel::MixtureAveraged ) {
    // attach requested components
    std::unique_ptr<IChemistry> chemistry( nullptr );
    switch( type ) {
    case Type::Basic :
      chemistry = std::unique_ptr<IChemistry>( new Thermo( mechanism ) );
      break;
    case Type::Transport:
      chemistry = std::unique_ptr<IChemistry>
        ( new Transport( new Thermo( mechanism )
                         , transport_model ) );   
      break;
    case Type::Kinetics:
      chemistry = std::unique_ptr<IChemistry>
        ( new Kinetics( new Thermo( mechanism ) ) );
      break;
    case Type::TransportAndKinetics:
      chemistry = std::unique_ptr<IChemistry>
        ( new Kinetics( new Transport( new Thermo(mechanism)
                                       , transport_model ) ) );
      break;
    default:
      return nullptr;
    }

    assert( chemistry );
    
    // create chemistry
    chemistry->Create();

    // can simply return chemistry instead, move ctor would be used
    return std::move(chemistry);
  }

} // namespace Chemistry

#endif //  CHEMISTRY_HPP
