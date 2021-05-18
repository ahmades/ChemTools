#include <vector>
#include <limits>

#include <cmath>

#include <catch2/catch.hpp>

#include "config.hpp"
#include "test_config.hpp"
#include "chemistry.hpp"

struct Input {
  std::string mechanism;
  Chemistry::Type chemistry_type;
  Chemistry::TransportModel transport_model;
  double temperature;
  double pressure;
  Cantera::compositionMap composition;

  Input( std::string const mechanism_
         , Chemistry::Type const chemistry_type_
         , Chemistry::TransportModel const transport_model_
         , double const temperature_
         , double pressure_
         , Cantera::compositionMap const& composition_ )
    : mechanism( mechanism_ )
    , chemistry_type{ chemistry_type_ }
    , transport_model{ transport_model_ }
    , temperature{ temperature_ }
    , pressure{ pressure_ }
    , composition( composition_ )
  {}
  
  void SetChemistryType( Chemistry::Type const chemistry_type_ ) {
    chemistry_type = chemistry_type_;
  }

  void SetTransportModel( Chemistry::TransportModel const transport_model_ ) {
    transport_model = transport_model_;
  }
};

struct ThermoOutput {
  double temperature;
  double pressure;
  double density;
  double molar_enthalpy;
  double molar_entropy;
  double molar_specific_heat;
};

struct TransOutput {
  double viscosity;
  double thermal_conductivity;
};

struct RateOfProgress {
  std::vector<double> fwd;
  std::vector<double> rev;
  std::vector<double> net;
  RateOfProgress()
    : fwd()
    , rev()
    , net()
  {}
  RateOfProgress( size_t n )
    : fwd( n )
    , rev( n )
    , net( n )
  {}
};

struct KinOutput {
  RateOfProgress rop;
};

static std::unique_ptr<Chemistry::IChemistry> Equilibrium( Input const& input ) {
 
  std::unique_ptr<Chemistry::IChemistry> chemistry{ nullptr };
  try {
    
    // instantiate chemistry object
    chemistry = Create( input.mechanism
                        , input.chemistry_type
                        , input.transport_model );
    assert( chemistry );

    // get thermo object
    std::unique_ptr<Cantera::ThermoPhase>& thermo = chemistry->GetThermo();
    assert( thermo );

    // set thermo state 
    thermo->setState_TPX( input.temperature
                          , input.pressure
                          , input.composition );

    // compute equilibrium solution
    thermo->equilibrate( "HP" );
    
  } catch( Cantera::CanteraError const& err ) {
    std::cerr << err.what() << std::endl;
  }

  return std::move( chemistry );
}

static ThermoOutput TestThermo( std::unique_ptr<Cantera::ThermoPhase> const& thermo
                                , bool verbose = false ) {
  ThermoOutput output;
  output.temperature         = thermo->temperature();
  output.pressure            = thermo->pressure();
  output.density             = thermo->density();
  output.molar_enthalpy      = thermo->enthalpy_mole();
  output.molar_entropy       = thermo->entropy_mole();
  output.molar_specific_heat = thermo->cp_mole();

  if( verbose ) {
    Cantera::writelog( "\nThermo properties at the equilibrium state:\n"
                       "Temperature:          {:14.5g} K\n"
                       "Pressure:             {:14.5g} Pa\n"
                       "Density:              {:14.5g} kg/m3\n"
                       "Molar Enthalpy:       {:14.5g} J/kmol\n"
                       "Molar Entropy:        {:14.5g} J/kmol-K\n"
                       "Molar cp:             {:14.5g} J/kmol-K\n"
                       , output.temperature
                       , output.pressure
                       , output.density
                       , output.molar_enthalpy
                       , output.molar_entropy
                       , output.molar_specific_heat );
  }
  
  return output;
}

static TransOutput TestTrans( std::unique_ptr<Cantera::Transport> const& trans
                              , bool verbose = false  ) {  
  TransOutput output;
  output.viscosity            = trans->viscosity();
  output.thermal_conductivity = trans->thermalConductivity();

  if( verbose ) {
    Cantera::writelog( "\nTransport properties at the equilibrium state:\n"
                       "Viscosity:            {:14.5g} kg/m-s\n"
                       "Thermal Conductivity: {:14.5g} W/m-K\n"
                       , output.viscosity
                       , output.thermal_conductivity );
  }
  
  return output;
}

static KinOutput TestKin( std::unique_ptr<Cantera::Kinetics> const& kinetics
                          , bool verbose = false  ) { 
  KinOutput output;
  size_t const n_rxns = kinetics->nReactions();
  
  output.rop = RateOfProgress( n_rxns );
  kinetics->getFwdRatesOfProgress( output.rop.fwd.data() );
  kinetics->getRevRatesOfProgress( output.rop.rev.data() );
  kinetics->getNetRatesOfProgress( output.rop.net.data() );

  if( verbose ) {
    Cantera::writelog("\nReactions and their forward, reverse and net rates of progress:\n");
    for( size_t i = 0; i < n_rxns; i++ ) {
      Cantera::writelog( "{:6s} {:35s} {:14.5g} {:14.5g} {:14.5g}  kmol/m3/s\n"
                         , 'R' + std::to_string( i + 1 )
                         , kinetics->reactionString( i )
                         , output.rop.fwd[i]
                         , output.rop.rev[i]
                         , output.rop.net[i] );
    }
  }

  return output;
}

static void TestThermoOutput( ThermoOutput const& output ) {
  CHECK( Approx(  1.706732e+03 ).epsilon( TestConfig::tolerance ) == output.temperature );
  CHECK( Approx(  2.026500e+05 ).epsilon( TestConfig::tolerance ) == output.pressure );
  CHECK( Approx(  3.239982e-01 ).epsilon( TestConfig::tolerance ) == output.density );
  CHECK( Approx( -5.626613e+06 ).epsilon( TestConfig::tolerance ) == output.molar_enthalpy );
  CHECK( Approx(  2.433695e+05 ).epsilon( TestConfig::tolerance ) == output.molar_entropy );
  CHECK( Approx(  3.734765e+04 ).epsilon( TestConfig::tolerance ) == output.molar_specific_heat );
}

static void TestTransOutput( TransOutput const& output ) {
  CHECK( Approx( 5.837424e-05 ).epsilon( TestConfig::tolerance ) == output.viscosity );
  CHECK( Approx( 1.747307e-01 ).epsilon( TestConfig::tolerance ) == output.thermal_conductivity );
}

static void TestKinOutput( KinOutput const& output ) {
  RateOfProgress const& rop = output.rop;
  for( size_t i = 0; i < rop.net.size(); i++ ) {
    CHECK_THAT( output.rop.net[i], Catch::Matchers::WithinAbs( 0.0, TestConfig::tolerance ) );
  }
}

SCENARIO( "Equilibtium state can be computed", "[chemistry]" ) {
 
  GIVEN( "Chemistry parameters without a chemistry type" ) {
    // input parameters
    std::string const mechanism = "../data/gri30.cti";
    Chemistry::Type const chemistry_type = Chemistry::Type::Invalid;
    Chemistry::TransportModel const transport_model = Chemistry::TransportModel::None;
    double const temperature = 500.0;
    double const pressure =  2.0 * Cantera::OneAtm;
    Cantera::compositionMap const composition = {
      { "CH4" , 1.00 },
      { "O2"  , 1.00 },
      { "N2"  , 3.76 }
    };

    // set input with chemistry_type set to invalid
    Input input( mechanism
                 , chemistry_type
                 , transport_model
                 , temperature
                 , pressure
                 , composition );

    try {
      
      WHEN( "Basic chemsitry is requested" ) {
        // reset chemistry type
        input.SetChemistryType( Chemistry::Type::Basic );
        // equilibrate
        std::unique_ptr<Chemistry::IChemistry> const chemistry( Equilibrium( input ) );
      
        // test output
        THEN("Check thermo properties") {
          // get thermo object
          std::unique_ptr<Cantera::ThermoPhase> const& thermo = chemistry->GetThermo();
          assert( thermo );
          TestThermoOutput( TestThermo( thermo, TestConfig::verbose ) );
        }
      }

      WHEN( "Basic chemsitry with transport is requested" ) {
        // reset chemistry type
        input.SetChemistryType( Chemistry::Type::Transport );
        input.SetTransportModel( Chemistry::TransportModel::MixtureAveraged );
      
        // equilibrate
        std::unique_ptr<Chemistry::IChemistry> const chemistry( Equilibrium( input ) );
      
        // test output
        THEN("Check thermo and transport properties") {
          {
            // get thermo object
            std::unique_ptr<Cantera::ThermoPhase> const& thermo = chemistry->GetThermo();
            //assert( thermo );
            REQUIRE_FALSE( thermo == nullptr );
            TestThermoOutput( TestThermo( thermo, TestConfig::verbose ) );
          }

          {
            // get transport object
            std::unique_ptr<Cantera::Transport> const& trans = chemistry->GetTrans();
            REQUIRE_FALSE( trans == nullptr );
            TestTransOutput( TestTrans( trans, TestConfig::verbose ) );
          }
        }
      }

      WHEN( "Basic chemsitry with kinetics and kinetics is requested" ) {
        // reset chemistry type
        input.SetChemistryType( Chemistry::Type::Kinetics );
      
        // equilibrate
        std::unique_ptr<Chemistry::IChemistry> const chemistry( Equilibrium( input ) );
      
        // test output
        THEN("Check thermo and kinetic properties") {
          {
            // get thermo object
            std::unique_ptr<Cantera::ThermoPhase> const& thermo = chemistry->GetThermo();
            REQUIRE_FALSE( thermo == nullptr );
            TestThermoOutput( TestThermo( thermo, TestConfig::verbose ) );
          }

          {
            // get kinetics object
            std::unique_ptr<Cantera::Kinetics> const& kinetics = chemistry->GetKinetics();
            REQUIRE_FALSE( kinetics == nullptr );
            TestKinOutput( TestKin( kinetics, TestConfig::verbose ) );
          }
        }
      
      }
    
      WHEN( "Basic chemsitry with transport and kinetics is requested" ) {
        // reset chemistry type
        input.SetChemistryType( Chemistry::Type::TransportAndKinetics );
        input.SetTransportModel( Chemistry::TransportModel::MixtureAveraged );
      
        // equilibrate
        std::unique_ptr<Chemistry::IChemistry> const chemistry( Equilibrium( input ) );
      
        // test output
        THEN("Check thermo, transport and kinetic properties") {
          {
            // get thermo object
            std::unique_ptr<Cantera::ThermoPhase> const& thermo = chemistry->GetThermo();
            REQUIRE_FALSE( thermo == nullptr );
            TestThermoOutput( TestThermo( thermo, TestConfig::verbose ) );
          }

          {
            // get transport object
            std::unique_ptr<Cantera::Transport> const& trans = chemistry->GetTrans();
            REQUIRE_FALSE( trans == nullptr );
            TestTransOutput( TestTrans( trans, TestConfig::verbose ) );
          }

          {
            // get kinetics object
            std::unique_ptr<Cantera::Kinetics> const& kinetics = chemistry->GetKinetics();
            REQUIRE_FALSE( kinetics == nullptr );
            TestKinOutput( TestKin( kinetics, TestConfig::verbose ) );
          }
        }
      
      }

    } catch( Cantera::CanteraError const& err ) {
      std::cout << err.what() << std::endl;
    }

    Chemistry::CleanUp();
    
  }
}


