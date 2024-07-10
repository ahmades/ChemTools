#ifndef APP_REACTORS_BASE_HPP
#define APP_REACTORS_BASE_HPP

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "chemistry.hpp"
#include "embedded_python/functions.hpp"
#include "results/observer.hpp"
#include "sundials/cvode/interface.hpp"

namespace Apps
{

  namespace Reactor
  {

    enum class Type
    {
      Const_Pres,
      Const_Vol,
      Const_Temp_Pres,
      Const_Temp_Vol
    };

    // ------------------
    // --- class Base ---
    // ------------------

    class Base : public SUNDIALS::CVODE::Client, public Results::Subject
    {

    public:
      Base();

      Base(Cantera::ThermoPhase& thermo,
           Cantera::Kinetics& kinetics,
           double const temperature,
           double const pressure,
           std::vector<double> const& mass_fractions,
           double const relative_solver_tolerance,
           double const absolute_solver_tolerance,
           double const total_simulation_time,
           bool const write_results);

      virtual ~Base() = default;

      Base(Base const&) = delete;

      Base& operator=(Base const&) = delete;

      double Temperature() const;

      double Pressure() const;

      double Density() const;

      std::vector<double> const& MassFractions() const;

      int RightHandSide(realtype const /*time*/,
                        realtype* const state,
                        realtype* rhs) override;

      void StoreResults(double const time);

    protected:
      Cantera::ThermoPhase& m_thermo;
      Cantera::Kinetics& m_kinetics;
      size_t m_nspecs;
      double m_density;                           // unit = [kg/m^3]
      double m_temperature;                       // unit = [K]
      double m_pressure;                          // unit = [Pa] = [kg/m/s^2]
      std::vector<double> m_mass_fractions;       // unit = [-]
      std::vector<double> m_net_production_rates; // unit = [kmol/m^3/s]

      void SetInitialState();

      virtual void SetThermoState(realtype const* const state) = 0;

      void SpeciesDerivative(realtype* const derivative);

    private:
      std::vector<double> m_molecular_weights; // unit = [kg/kmol]
      double m_time;
      double m_relative_solver_tolerance;
      double m_absolute_solver_tolerance;
      double m_total_simulation_time;
      bool m_write_results;

      void CompleteInstantiation();

      void SetSolverParameters();

      void RegisterBasicResults();

      virtual void UpdateSolution();
    };

    // ------------------------------------------------------------------------
    // --- class EnergyEnabled
    // ------------------------------------------------------------------------

    class EnergyEnabled : public Base
    {
    public:
      EnergyEnabled(Cantera::ThermoPhase& thermo,
                    Cantera::Kinetics& kinetics,
                    double const temperature,
                    double const pressure,
                    std::vector<double> const& mass_fractions,
                    double const relative_solver_tolerance,
                    double const absolute_solver_tolerance,
                    double const total_simulation_time,
                    bool const write_results);

      virtual ~EnergyEnabled() = default;

      EnergyEnabled(EnergyEnabled const&) = delete;

      EnergyEnabled& operator=(EnergyEnabled const&) = delete;

      void UpdateSolution() override;

      int RightHandSide(realtype const /*time*/,
                        realtype* const state,
                        realtype* rhs) override;

    protected:
      std::vector<double> m_energy; // enthalpy or internal energy [J/kmol]

      void SetInitialState();

      virtual void SetThermoState(realtype const* const state) = 0;

      virtual double EnergyDerivative() = 0;

    private:
      void RegisterAdditionalResults();
    };

  } // namespace Reactor

} // namespace Apps

#endif // APP_REACTORS_BASE_HPP
