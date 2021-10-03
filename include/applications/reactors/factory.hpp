#ifndef APP_REACTORS_FACTORY_HPP
#define APP_REACTORS_FACTORY_HPP

#include <memory>

#include "applications/reactors/base.hpp"
#include "applications/reactors/reactors.hpp"

namespace Apps {
  
  namespace Reactor {

    inline std::unique_ptr<Base> Create( Type const reactor_type
                                         , Cantera::ThermoPhase& thermo
                                         , Cantera::Kinetics& kinetics
                                         , double const temperature
                                         , double const pressure
                                         , std::vector<double> const& mass_fractions
                                         , double const relative_solver_tolerance
                                         , double const absolute_solver_tolerance
                                         , double const total_simulation_time
                                         , bool const write_results ) {
      switch( reactor_type ) {
      case Type::Const_Pres :
        return std::make_unique<ConstPres>( thermo
                                            , kinetics
                                            , temperature
                                            , pressure
                                            , mass_fractions
                                            , relative_solver_tolerance
                                            , absolute_solver_tolerance
                                            , total_simulation_time
                                            , write_results);
      case Type::Const_Vol :
        return std::make_unique<ConstVol>( thermo
                                           , kinetics
                                           , temperature
                                           , pressure
                                           , mass_fractions
                                           , relative_solver_tolerance
                                           , absolute_solver_tolerance
                                           , total_simulation_time
                                           , write_results );
      case Type::Const_Temp_Pres :
        return std::make_unique<ConstTempPres>( thermo
                                                , kinetics
                                                , temperature
                                                , pressure
                                                , mass_fractions
                                                , relative_solver_tolerance
                                                , absolute_solver_tolerance
                                                , total_simulation_time
                                                , write_results );
      case Type::Const_Temp_Vol :
        return std::make_unique<ConstTempVol>( thermo
                                               , kinetics
                                               , temperature
                                               , pressure
                                               , mass_fractions
                                               , relative_solver_tolerance
                                               , absolute_solver_tolerance
                                               , total_simulation_time
                                               , write_results );
      default :
        throw( "Invalid reactor" );
        return nullptr;
      }
    }
    
  } // namespace Reactor
  
} // namespace Apps

#endif // APP_REACTORS_FACTORY_HPP
