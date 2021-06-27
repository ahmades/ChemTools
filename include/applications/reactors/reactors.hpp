#ifndef APP_REACTORS_CONST_PROP_HPP
#define APP_REACTORS_CONST_PROP_HPP

#include "applications/reactors/base.hpp"

namespace Apps {
  
  namespace Reactor {

    // --- class ConstPres
    
    class ConstPres : public EnergyEnabled {
    public:
      ConstPres( Cantera::ThermoPhase* const thermo
                 , Cantera::Kinetics* const kinetics
                 , double const temperature
                 , double const pressure
                 , std::vector<double> const& mass_fractions
                 , double const relative_solver_tolerance
                 , double const absolute_solver_tolerance
                 , double const total_simulation_time
                 , bool const write_results );
      ~ConstPres() {}
    private:
      void SetThermoState( realtype const * const state ) override;
      double EnergyDerivative() override;
    };

    // --- class ConstVol
    
    class ConstVol : public EnergyEnabled {
    public:
      ConstVol( Cantera::ThermoPhase* const thermo
                , Cantera::Kinetics* const kinetics
                , double const temperature
                , double const pressure
                , std::vector<double> const& mass_fractions
                , double const relative_solver_tolerance
                , double const absolute_solver_tolerance
                , double const total_simulation_time
                , bool const write_results );
      ~ConstVol() {}
    private:
      void SetThermoState( realtype const * const state ) override;
      double EnergyDerivative() override;
    };

    // --- class ConstTempPres
    
    class ConstTempPres : public Base {
    public:
      ConstTempPres( Cantera::ThermoPhase* const thermo
                     , Cantera::Kinetics* const kinetics
                     , double const temperature
                     , double const pressure
                     , std::vector<double> const& mass_fractions
                     , double const relative_solver_tolerance
                     , double const absolute_solver_tolerance
                     , double const total_simulation_time
                     , bool const write_results );
      ~ConstTempPres() {}
    private:
      void SetThermoState( realtype const * const state ) override;
    };

    // --- class ConstTempVol
    
    class ConstTempVol : public Base {
    public:
      ConstTempVol( Cantera::ThermoPhase* const thermo
                    , Cantera::Kinetics* const kinetics
                    , double const temperature
                    , double const pressure
                    , std::vector<double> const& mass_fractions
                    , double const relative_solver_tolerance
                    , double const absolute_solver_tolerance
                    , double const total_simulation_time
                    , bool const write_results );
      ~ConstTempVol() {}
    private:
      void SetThermoState( realtype const * const state ) override;
    };

  } // namespace Reactor
  
} // namespace Apps

#endif // APP_REACTORS_CONST_PROP_HPP
