#include <limits>
#include <vector>

#include "applications/reactors/base.hpp"

namespace Apps
{

  namespace Reactor
  {

    // ------------------------------------------------------------------------
    // --- class Base
    // ------------------------------------------------------------------------

    Base::Base(Cantera::ThermoPhase& thermo, Cantera::Kinetics& kinetics, double const temperature, double const pressure, std::vector<double> const& mass_fractions, double const relative_solver_tolerance, double const absolute_solver_tolerance, double const total_simulation_time, bool const write_results)
        : m_thermo(thermo), m_kinetics(kinetics), m_nspecs(m_thermo.nSpecies()), m_density{std::numeric_limits<double>::signaling_NaN()}, m_temperature{temperature}, m_pressure{pressure}, m_mass_fractions(mass_fractions), m_net_production_rates(m_nspecs, std::numeric_limits<double>::signaling_NaN()), m_molecular_weights(m_nspecs, std::numeric_limits<double>::signaling_NaN()), m_time{std::numeric_limits<double>::signaling_NaN()}, m_relative_solver_tolerance{relative_solver_tolerance}, m_absolute_solver_tolerance{absolute_solver_tolerance}, m_total_simulation_time{total_simulation_time}, m_write_results{write_results}
    {
      CompleteInstantiation();
    }

    double Base::Temperature() const
    {
      return m_temperature;
    }

    double Base::Pressure() const
    {
      return m_pressure;
    }

    double Base::Density() const
    {
      return m_density;
    }

    std::vector<double> const& Base::MassFractions() const
    {
      return m_mass_fractions;
    }

    int Base::RightHandSide(realtype const /*time*/
                            ,
                            realtype* const state,
                            realtype* rhs)
    {
      SetThermoState(state);
      SpeciesDerivative(rhs);
      return 0;
    }

    void Base::StoreResults(double const time)
    {
      // update time
      m_time = time;
      // update local state
      UpdateSolution();
      // notify the observer
      if (m_write_results)
        {
          Results::Subject::NotifyObserver();
        }
    }

    void Base::SetInitialState()
    {
      SUNDIALS::CVODE::Client::SetState(m_mass_fractions);
    }

    void Base::SpeciesDerivative(realtype* const derivative)
    {
      m_kinetics.getNetProductionRates(m_net_production_rates.data());
      for (size_t i = 0; i < m_nspecs; ++i)
        {
          // [kmol/(m^3.s)] * [kg/kmol] / [kg/m3] = [1/s] = unit_of( dY/dt )
          derivative[i] = m_net_production_rates[i]
                          * m_molecular_weights[i] / m_density;
        }
    }

    void Base::CompleteInstantiation()
    {
      // compute molecular weights
      m_thermo.getMolecularWeights(m_molecular_weights.data());

      // initialise the client's simulation parameters
      SetSolverParameters();

      // register time and mass farctions results
      if (m_write_results)
        {
          RegisterBasicResults();
        }
    }

    void Base::SetSolverParameters()
    {
      // set the tolerances
      SUNDIALS::CVODE::Client::SetTolerances(m_relative_solver_tolerance, m_absolute_solver_tolerance);
      // set the integration start and stop times
      SUNDIALS::CVODE::Client::SetIntegrationTime(0, m_total_simulation_time);
    }

    void Base::RegisterBasicResults()
    {

      // TODO: result file name must be constructed from case parameters
      Results::Subject::results_meta.SetFileName("test");

      Results::Subject::results_meta.Register("Time", "Time", "Time", "t", "s", &m_time);

      for (size_t i = 0; i < m_nspecs; ++i)
        {
          std::string const species_name = m_thermo.speciesName(i);
          std::string const notation = "Y_" + m_thermo.speciesName(i);
          Results::Subject::results_meta.Register("Mass fractions", species_name, notation, notation, "-", &m_mass_fractions[i]);
        }
    }

    void Base::UpdateSolution()
    {
      double const* const state_ptr = Client::State();
      std::copy(state_ptr, state_ptr + m_nspecs, m_mass_fractions.begin());
    }

    // ------------------------------------------------------------------------
    // --- class EnergyEnabled
    // ------------------------------------------------------------------------

    EnergyEnabled::EnergyEnabled(Cantera::ThermoPhase& thermo, Cantera::Kinetics& kinetics, double const temperature, double const pressure, std::vector<double> const& mass_fractions, double const relative_solver_tolerance, double const absolute_solver_tolerance, double const total_simulation_time, bool const write_results)
        : Base(thermo, kinetics, temperature, pressure, mass_fractions, relative_solver_tolerance, absolute_solver_tolerance, total_simulation_time, write_results), m_energy(m_nspecs, std::numeric_limits<double>::signaling_NaN())
    {
      // register additional results (temperature)
      if (write_results)
        {
          RegisterAdditionalResults();
        }
    }

    void EnergyEnabled::UpdateSolution()
    {
      double const* const state_ptr = Client::State();
      m_temperature = state_ptr[0];
      std::copy(state_ptr + 1, state_ptr + 1 + m_nspecs, m_mass_fractions.begin());
    }

    int EnergyEnabled::RightHandSide(realtype const /*time*/
                                     ,
                                     realtype* const state,
                                     realtype* rhs)
    {
      SetThermoState(state);
      SpeciesDerivative(rhs + 1);
      rhs[0] = EnergyDerivative();
      return 0;
    }

    void EnergyEnabled::SetInitialState()
    {
      // set the initial state
      std::vector<double> initial_state(m_mass_fractions);
      initial_state.insert(initial_state.begin(), m_temperature);
      SUNDIALS::CVODE::Client::SetState(initial_state);
    }

    void EnergyEnabled::RegisterAdditionalResults()
    {
      Results::Subject::results_meta.Register("Temperature", "Temperature", "Temperature", "T", "K", &m_temperature);
    }

  } // namespace Reactor

} // namespace Apps
