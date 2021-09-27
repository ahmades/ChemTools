#include <vector>
#include <limits>

#include <catch2/catch.hpp>

#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/color.h>

#include "config.hpp"
#include "test_config.hpp"
#include "chemistry.hpp"
#include "results/hdf5_writer.hpp"
#include "applications/reactors/factory.hpp"

using namespace Apps;

TEST_CASE( "Test reactors", "[reactors]" ) {

  std::string const mechanism = "../data/gri30.cti";
  std::unique_ptr<Chemistry::IChemistry> chemistry =
    Chemistry::Create( mechanism, Chemistry::Type::Kinetics );
  Cantera::ThermoPhase* const thermo = chemistry->ThermoPtr();
  Cantera::Kinetics* const kinetics = chemistry->KineticsPtr();

  // target: http://combustion.berkeley.edu/gri-mech/version30/targets30/ig.1b.html
  /*double pressure = 2.04 * Cantera::OneAtm;
  double temperature = 1700.0;
  std::vector<double> mole_fractions( thermo->nSpecies(), 0.0 );
  mole_fractions[thermo->speciesIndex("CH4")] = 0.091;
  mole_fractions[thermo->speciesIndex("O2")]  = 0.182;
  mole_fractions[thermo->speciesIndex("Ar")]  = 0.727;
  thermo->setState_TPX( temperature, pressure, mole_fractions.data()  );
  std::vector<double> mass_fractions( thermo->nSpecies(), 0.0 );
  thermo->getMassFractions( mass_fractions.data() );*/
  
  // target: http://combustion.berkeley.edu/gri-mech/version30/targets30/ig.st1a.html
  /*double pressure = 6.1 * Cantera::OneAtm;
  double temperature = 1356.0;
  std::vector<double> mole_fractions( thermo->nSpecies(), 0.0 );
  mole_fractions[thermo->speciesIndex("CH4")]  = 0.0329;
  mole_fractions[thermo->speciesIndex("C2H6")] = 0.0021;
  mole_fractions[thermo->speciesIndex("O2")]   = 0.0700;
  mole_fractions[thermo->speciesIndex("Ar")]   = 0.8950;
  thermo->setState_TPX( temperature, pressure, mole_fractions.data()  );
  std::vector<double> mass_fractions( thermo->nSpecies(), 0.0 );
  thermo->getMassFractions( mass_fractions.data() );*/

  // cantera
  double const pressure = Cantera::OneAtm;
  double const temperature = 1001.0;
  size_t const n_specs = thermo->nSpecies();
  std::vector<double> mole_fractions( n_specs, 0.0 );
  mole_fractions[thermo->speciesIndex("H2")] = 2.0;
  mole_fractions[thermo->speciesIndex("O2")] = 1.0;
  mole_fractions[thermo->speciesIndex("N2")] = 4.0;
  thermo->setState_TPX( temperature, pressure, mole_fractions.data()  );
  std::vector<double> mass_fractions( n_specs, 0.0 );
  thermo->getMassFractions( mass_fractions.data() );

  double const total_simulation_time = 1.0e-3;
  double const relative_solver_tolerance = 1.0e-9;
  double const absolute_solver_tolerance = 1.0e-15;
  bool const write_results = true;
  
  std::unique_ptr<Reactor::Base> reactor =
    Reactor::Create( Reactor::Type::Const_Vol
                     , thermo
                     , kinetics
                     , temperature
                     , pressure
                     , mass_fractions
                     , relative_solver_tolerance
                     , absolute_solver_tolerance
                     , total_simulation_time
                     , write_results );

  // attach to result observer
  if( write_results ) {
    Results::HDF5Writer result_writer( *(reactor.get()), "." );
  }

  
  // select solution strategy
  SUNDIALS::CVODE::Dense strategy; 
  
  // set strategy
  SUNDIALS::CVODE::Solver cvode( &strategy );

  // set linear multistep method
  cvode.SetLinearMultiStepMethod( SUNDIALS::CVODE::Types::LinearMultisptepMethod::BDF );
  
  // set client data
  cvode.SetClientData( reactor.get() );

  // set additional solver options
  cvode.SetMaxNumSteps( 20000 );
  
  // initialise
  REQUIRE( 0 == cvode.Initialise() );

 // integrate
  {
    realtype time = 0.0;
    realtype const time_step = 1.0e-5; //1.0e-7;
    bool ignited = false;
    realtype time_ignition = 0.0;
    
    while( time < total_simulation_time ) {
      // set target time
      time += time_step;

      // integrate
      int const ret_code = cvode.Integrate( time );
      if( ret_code < CV_SUCCESS ) {
        fmt::print( "Solver returned {}.\n", ret_code );
      }
      
      if( time < total_simulation_time ) {
        REQUIRE( CV_SUCCESS == ret_code );
      } else {
        REQUIRE( CV_TSTOP_RETURN == ret_code );
      }
      
      reactor->StoreResults( time );
      
      double temp = cvode.Solution()[0];
      double mf = cvode.Solution()[1 + thermo->speciesIndex("H2")];
      if( temp >= 1756.0 && ignited == false ) {
        time_ignition = time;
        ignited = true;
      }
      
      if( TestConfig::verbose ) {
        fmt::print( "t = {} s T = {} K Y_H2 = {} \n"
                    , time
                    , temp
                    , mf );
      }
    }

    if( TestConfig::verbose ) {
      fmt::print( "\nIgnition at t = {} ms.\n",  time_ignition * 1.0e+3 );
      
      std::string const title = fmt::format( fmt::emphasis::underline
                                             | fmt::emphasis::bold
                                             , "Final statistics:" );    
      fmt::print( "{}\n{}\n"
                  , title
                  , cvode.PrintSolverStatistics() );
    }
  }
}
  
