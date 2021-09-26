#ifndef CHEMISTRY_HPP
#define CHEMISTRY_HPP

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
  void CleanUp();

  //
  // -------------- chemisty interface
  //
  class IChemistry {
    friend class Decorator;
  public:
    IChemistry();
    virtual ~IChemistry() = default;
    virtual void Create() = 0;
    virtual Cantera::ThermoPhase* ThermoPtr();
    virtual Cantera::Transport* TransPtr();
    virtual Cantera::Kinetics* KineticsPtr();
  protected:
    std::unique_ptr<Cantera::ThermoPhase> thermo;
    std::unique_ptr<Cantera::Transport> trans;
    std::unique_ptr<Cantera::Kinetics> kinetics;
    virtual std::unique_ptr<Cantera::Transport>& TransUPtr();
    virtual std::unique_ptr<Cantera::Kinetics>& KineticsUPtr();
  };

  //
  // -------------- simplest IChemistry instance only has thermo
  //
  class Thermo: public IChemistry {
  public:
    explicit Thermo( std::string const& mechanism_ );
    ~Thermo() = default;
    void Create() override;
  private:
    std::string mechanism;
  };

  //
  // -------------- decorator
  //
  class Decorator: public IChemistry {
  public:
    explicit Decorator( std::unique_ptr<IChemistry> chemistry_ );
    ~Decorator() = default; 
    virtual void Create();
    Cantera::ThermoPhase* ThermoPtr() override;
    Cantera::Transport* TransPtr() override;
    Cantera::Kinetics* KineticsPtr() override;
  protected:
    std::unique_ptr<Cantera::Transport>& TransUPtr() override;
    std::unique_ptr<Cantera::Kinetics>& KineticsUPtr() override;
  private:
    std::unique_ptr<IChemistry> chemistry;
  };
  
  //
  // -------------- attaches transport
  //
  class Transport : public Decorator {
  public:
    explicit Transport( std::unique_ptr<IChemistry> chemistry
                        , TransportModel model_ = TransportModel::MixtureAveraged );
    ~Transport() = default;
    void Create() override;
  private:
    TransportModel model;
  };

  //
  // -------------- attaches kinetics
  //
  class Kinetics : public Decorator {
  public:
    explicit Kinetics( std::unique_ptr<IChemistry> chemistry );
    ~Kinetics() = default;
    void Create() override;
  };

  //
  // -------------- convenience factory
  //
  std::unique_ptr<IChemistry> Create( std::string const& mechanism
                                      , Type type
                                      , TransportModel transport_model
                                      = TransportModel::MixtureAveraged );

} // namespace Chemistry

#endif //  CHEMISTRY_HPP
