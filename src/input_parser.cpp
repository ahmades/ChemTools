#include <string>
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_set>

#include <cassert>
#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>

#include "yaml-cpp/yaml.h"

#include <fmt/format.h>

#include "input/types.hpp"
#include "input/parser.hpp"
#include "chemistry.hpp"

namespace YAML {

  // decoder of case sets
  template<>
  struct convert<Input::CaseSet> {
    static bool decode( Node const& node
			, Input::CaseSet& case_set ) {
      
      assert( node.IsMap() );

      // mechanism
      {
	Node const& node_mechanism = node["Mechanism"];
	assert( node_mechanism.IsMap());
	Input::Mechanism mechanism;
	mechanism.SetPath( node_mechanism["Path"].as<std::string>() );
	if( node_mechanism["Description"].IsDefined() ) {
	  mechanism.SetDescription( node_mechanism["Description"].as<std::string>() );
	}
	case_set.SetMechanism( mechanism );
      }

      // description
      {
	Node const& node_description = node["Description"];
        if( node_description.IsDefined() ) {
          assert( node_description.IsScalar());
	  case_set.SetDescription( node_description.as<std::string>() );
	}
      }
      
      // cases
      {
	Node const& node_cases = node["Cases"];
	assert( node_cases.IsMap() );
	for( YAML::const_iterator it = node_cases.begin(); it != node_cases.end(); ++it) {
	  YAML::Node const& node_key = it->first;
	  YAML::Node const& node_value = it->second;
	  if( node_key.as<std::string>() == "Case" && node_value.IsMap() ) {
	    Input::Case cas;

            // description
	    {
	      YAML::Node const& node = node_value["Description"];
              if( node.IsDefined() ) {
                assert( node.IsScalar() );
                cas.SetDescription( node.as<std::string>() );
              }
	    }

            // pressure
	    {
	      YAML::Node const& node = node_value["Pressure"];
	      assert( node.IsMap() );
	      /*cas.SetPressure( node["Value"].as<double>()
                , node["Unit"].as<std::string>() );*/
              cas.pressure.Set( node["Value"].as<double>()
                                , node["Unit"].as<std::string>() );
	    }

            // temperature
	    {
	      YAML::Node const& node = node_value["Temperature"];
	      assert( node.IsMap() );
              cas.temperature.Set( node["Value"].as<double>()
                                   , node["Unit"].as<std::string>() );
	    }

            // composition
	    {
	      YAML::Node const& node = node_value["Composition"];
	      assert( node.IsMap() );
              std::string const type = node["Type"].as<std::string>();
              if( "concentration" == type ) {
                cas.composition.Set( node["Value"].as<Input::CompositionPairs>()
                                     , type
                                     , node["Unit"].as<std::string>() );
              } else {
                cas.composition.Set( node["Value"].as<Input::CompositionPairs>()
                                     , node["Type"].as<std::string>() );
              }
	    }

            // time
	    {
	      YAML::Node const& node = node_value["Time"];
	      assert( node.IsMap() );
              cas.time.Set( node["Value"].as<double>()
                            , node["Unit"].as<std::string>() );
	    }
            
	    case_set.AddCase( cas );
	  }
	}
      }
      return true;
    }
  };
  
}

namespace Input {

  Parser::Parser()
    : config_file()
    , parameters()
  {}
  
  Parser::Parser( std::string const& config_file_ )
    : config_file( config_file_ )
    , parameters()
  {
    Parse();
    Validate();
  }

  Parameters const& Parser::GetParameters() const {
    return parameters;
  }
  
  void Parser::Parse() {
    Input::Reactor& reactor = parameters.reactor;
    YAML::Node const& root = YAML::LoadFile( config_file );
    for( YAML::Node const& node : root ) {
      if (node["Reactor"].as<std::string>() == "CONSTANT_PRESSURE") {
        reactor.const_P.push_back(node.as<CaseSet>());
      }
      if (node["Reactor"].as<std::string>() == "CONSTANT_VOLUME") {
        reactor.const_V.push_back(node.as<CaseSet>());
      }
      if (node["Reactor"].as<std::string>() == "CONSTANT_TEMPERATURE_PRESSURE") {
        reactor.const_TP.push_back(node.as<CaseSet>());
      }
      if (node["Reactor"].as<std::string>() == "CONSTANT_TEMPERATURE_VOLUME") {
        reactor.const_TV.push_back(node.as<CaseSet>());
      }
    }
  }

  static void MechanismPathIsValid( std::string const& inp_path_str ) {
    /* check:
       path is valid
       TODO: mechanism has at least thermo and kinetics, trans not needed 
     */
    boost::filesystem::path const path( inp_path_str );
    boost::system::error_code error_code;
    if( !boost::filesystem::exists( path, error_code ) ) {
      BOOST_FILESYSTEM_THROW( boost::filesystem::filesystem_error( "Input error"
                                                                   , path
                                                                   , error_code ) );
    }
  }

  static void TemperatureIsValid( Cantera::ThermoPhase const& thermo
                                  , RealScalar const& temperature_meta ) {
    double const temperature_min = thermo.minTemp();
    double const temperature_max = thermo.maxTemp();
    double const temperature_value = temperature_meta.Value();
    std::string const& temperature_unit = units::to_string( temperature_meta.Unit() );
    if( temperature_value < 0 ) {
      throw std::runtime_error( fmt::format( "Input error: the specified teperature is below absolute zero. "
                                             "T = {} {}.\n"
                                             , temperature_value
                                             , temperature_unit ) );
    } else if( temperature_value < temperature_min || temperature_value > temperature_max ) {
      // warning is needed here, not an error
      // should use logger here instead of cout
      fmt::print( "Input warning: the specified temperature"
                  " T = {} {} is outside the mechanism's applicable"
                  " thermodynamic range [{}, {}] {}.\n"
                  , temperature_value
                  , temperature_unit
                  , temperature_min
                  , temperature_max
                  , temperature_unit );
    }
  }

  static void PressureIsValid( RealScalar const& pressure_meta ) {
    double const pressure = pressure_meta.Value();
    if( pressure < 0.0 ) {
      throw std::runtime_error( fmt::format( "Input error: the specified pressure is negative. "
                                             "P = {} {}.\n"
                                             , pressure
                                             , units::to_string( pressure_meta.Unit() ) ) );
    }
  }

  static void CompositionIsValid( Cantera::ThermoPhase const& thermo
                                  , Composition const& composition_meta ) {
    
    CompositionPairs const& composition_pairs = composition_meta.Pairs();
    
    // check if species are unique
    {
      std::unordered_set<std::string> unique_species_set;
      std::vector<std::string> duplicate_species;
      bool duplicate_species_found = false;
      for( auto const& pair : composition_pairs ) {
        if( unique_species_set.find( pair.first ) != unique_species_set.end() ) {
          duplicate_species_found = true;
          duplicate_species.push_back( pair.first );
        } else {
          unique_species_set.insert( pair.first );
        }
      }
      
      if( duplicate_species_found ) {
        std::string str{};
        std::for_each( duplicate_species.begin(), duplicate_species.end()
                       , [&str]( auto const& species ) {
                         str += fmt::format("' '{}' '", species ); } );
        throw std::runtime_error( fmt::format( "Input error: found duplicate species: [{}].\n"
                                               , str ) );
      }
    }
    
    // check if species exist in mechanism
    {
      std::string invalid_species{};
      bool invalid_species_found = false;
      for( auto const& pair : composition_pairs ) {
        std::string const& species = pair.first;
        if( -1 == static_cast<int>( thermo.speciesIndex( species ) ) ) {
          invalid_species += fmt::format("' '{}' '", species );
          invalid_species_found = true;
        }
      }
      
      if( invalid_species_found ) {
        throw std::runtime_error( fmt::format( "Input error: found invalid species: [{}].\n"
                                               , invalid_species ) );
      }
    } // existance scope
    
    // check if species are positive
    {
      std::string negative_species{};
      bool negative_species_found = false;
      for( auto const& pair : composition_pairs ) {
        if( pair.second < 0 ) {
          negative_species += fmt::format("' '{}:{}' '", pair.first, pair.second );
          negative_species_found = true;
        }
      }
      if( negative_species_found ) {
        throw std::runtime_error( fmt::format( "Input error: found negative species: [{}].\n"
                                               , negative_species ) );
      }
    } // positivity scope
    
    // check if the specified mass or mole fractions sum up to unity
    {
      CompositionInputType const type = composition_meta.Type();
      if(    CompositionInputType::MassFraction == type
          || CompositionInputType::MoleFraction == type ) {
        double const sum =
          std::accumulate( composition_pairs.begin()
                           , composition_pairs.end()
                           , 0.0
                           , []( double const previous, std::pair<std::string, double> const& pair) {
                             return previous + pair.second; } );
        // TODO: tolerance?
        if( std::abs( 1.0 - sum ) > 1.0e-12 ) {
          throw std::runtime_error( fmt::format( "Input error: the {} fractions do not sum up to unity. The sum is {}.\n"
                                                 , ( type == CompositionInputType::MassFraction ) ? "mass" : "mole"
                                                 , sum ) );
        } 
      }
    } // conservativeness scopre
  }

  static void ValidateCaseSet( std::vector<CaseSet> const& case_sets ) {
    if( !case_sets.empty() ) {
      for( auto case_set : case_sets ) {
        // check mechanism
        std::string const& mechanism_path = case_set.GetMechanism().Path();
        //try {
        MechanismPathIsValid( mechanism_path );
        // } catch( boost::filesystem::filesystem_error const& fs_err ) {
        //   fmt::print( "{}\n", fs_err.what() );
        //   std::exit( EXIT_FAILURE );
        // }
        
        // create chemistry with thermo
        Chemistry::Thermo chemistry( mechanism_path );
        // get thermo ptr
        Cantera::ThermoPhase const& thermo = chemistry.thermo();
        
        // loop over all cases in case set
        std::vector<Case> const& cases = case_set.Cases();
        for( auto const& cas : cases ) {
          // check temperature
          TemperatureIsValid( thermo, cas.temperature );
          // check pressure
          PressureIsValid( cas.pressure );
          // check composition
          CompositionIsValid( thermo, cas.composition );       
        }
      } // case sets loop
    }
  }
    
  void Parser::Validate() const {
    ValidateCaseSet( parameters.reactor.const_P );
    ValidateCaseSet( parameters.reactor.const_V );
    ValidateCaseSet( parameters.reactor.const_TP );
    ValidateCaseSet( parameters.reactor.const_TV );
  } // Parser::Validate
  
} // namespace Input
