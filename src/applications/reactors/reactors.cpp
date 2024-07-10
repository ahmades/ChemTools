
#include "applications/reactors/reactors.hpp"

namespace Apps
{

  namespace Reactor
  {

    // --- class ConstPres

    ConstPres::ConstPres(Cantera::ThermoPhase& thermo,
                         Cantera::Kinetics& kinetics,
                         double const temperature,
                         double const pressure,
                         std::vector<double> const& mass_fractions,
                         double const relative_solver_tolerance,
                         double const absolute_solver_tolerance,
                         double const total_simulation_time,
                         bool const write_results)
        : EnergyEnabled(thermo,
                        kinetics,
                        temperature,
                        pressure,
                        mass_fractions,
                        relative_solver_tolerance,
                        absolute_solver_tolerance,
                        total_simulation_time,
                        write_results)
    {
      // set the client's initial state
      SetInitialState();
    }

    void ConstPres::SetThermoState(realtype const* const state)
    {
      // set the thrmo state
      m_thermo.setState_TPY(state[0], m_pressure, state + 1);

      // update the density
      m_density = m_thermo.density();
    }

    double ConstPres::EnergyDerivative()
    {
      // get partial molar enthalpies
      m_thermo.getPartialMolarEnthalpies(m_energy.data());

      // [J/kmol] * [kmol/m^3/s] = [J//m^3/s]
      double const inner_prod = std::inner_product(m_energy.begin(),
                                                   m_energy.end(),
                                                   m_net_production_rates.begin(),
                                                   0.0);

      // [J//m^3/s] / [kg/m^3] / [J/kg/K] = [K/s] = unit_of( dT/dt )
      return (-inner_prod / m_density / m_thermo.cp_mass());
    }

    // --- class ConsVol

    ConstVol::ConstVol(Cantera::ThermoPhase& thermo,
                       Cantera::Kinetics& kinetics,
                       double const temperature,
                       double const pressure,
                       std::vector<double> const& mass_fractions,
                       double const relative_solver_tolerance,
                       double const absolute_solver_tolerance,
                       double const total_simulation_time,
                       bool const write_results)
        : EnergyEnabled(thermo,
                        kinetics,
                        temperature,
                        pressure,
                        mass_fractions,
                        relative_solver_tolerance,
                        absolute_solver_tolerance,
                        total_simulation_time,
                        write_results)
    {
      // set the client's initial state
      SetInitialState();

      // set the initial thermo state and compute the initial density, then hold it constant
      m_thermo.setState_TPY(m_temperature, m_pressure, m_mass_fractions.data());
      m_density = thermo.density();
    }

    void ConstVol::SetThermoState(realtype const* const state)
    {
      // set the temperature, density and mass fractions
      m_thermo.setState_TRY(state[0], m_density, state + 1);
      // update the pressure
      m_pressure = m_thermo.pressure();
    }

    double ConstVol::EnergyDerivative()
    {
      // get the partial internal energies
      m_thermo.getPartialMolarIntEnergies(m_energy.data());

      // [J/kmol] * [kmol/m^3/s] = [J//m^3/s]
      double const inner_prod = std::inner_product(m_energy.begin(),
                                                   m_energy.end(),
                                                   m_net_production_rates.begin(),
                                                   0.0);

      // [J//m^3/s] / [kg/m^3] / [J/kg/K] = [K/s] = unit_of( dK/dt )
      return (-inner_prod / m_density / m_thermo.cv_mass());
    }

    // --- class ConstTempPres

    ConstTempPres::ConstTempPres(Cantera::ThermoPhase& thermo,
                                 Cantera::Kinetics& kinetics,
                                 double const temperature,
                                 double const pressure,
                                 std::vector<double> const& mass_fractions,
                                 double const relative_solver_tolerance,
                                 double const absolute_solver_tolerance,
                                 double const total_simulation_time,
                                 bool const write_results)
        : Base(thermo,
               kinetics,
               temperature,
               pressure,
               mass_fractions,
               relative_solver_tolerance,
               absolute_solver_tolerance,
               total_simulation_time,
               write_results)
    {
      // set the client's initial state
      SetInitialState();

      // set the initial temeperaturem pressure and mass fractions
      m_thermo.setState_TPY(m_temperature,
                            m_pressure,
                            m_mass_fractions.data());
      // compute the initial density and hold it constant
      m_pressure = m_thermo.pressure();
      // set the initial temperature and pressure and hold them constant
    }

    void ConstTempPres::SetThermoState(realtype const* const state)
    {
      // set the temperature, pressure and mass fractions
      m_thermo.setState_TPY(m_temperature,
                            m_pressure,
                            state);
      // update the density
      m_density = m_thermo.density();
    }

    // --- class ConstTempVol

    ConstTempVol::ConstTempVol(Cantera::ThermoPhase& thermo,
                               Cantera::Kinetics& kinetics,
                               double const temperature,
                               double const pressure,
                               std::vector<double> const& mass_fractions,
                               double const relative_solver_tolerance,
                               double const absolute_solver_tolerance,
                               double const total_simulation_time,
                               bool const write_results)
        : Base(thermo,
               kinetics,
               temperature,
               pressure,
               mass_fractions,
               relative_solver_tolerance,
               absolute_solver_tolerance,
               total_simulation_time,
               write_results)
    {
      // set the client's initial state
      SetInitialState();

      // set the initial temeperaturem pressure and mass fractions
      m_thermo.setState_TPY(m_temperature,
                            m_pressure,
                            m_mass_fractions.data());
      // compute the initial density and hold it constant
      m_density = m_thermo.density();
    }

    void ConstTempVol::SetThermoState(realtype const* const state)
    {
      // set the temperature, density and mass fractions
      m_thermo.setState_TRY(m_temperature,
                            m_density,
                            state);
      // update the pressure
      m_pressure = m_thermo.pressure();
    }

  } // namespace Reactor

} // namespace Apps
