#include <limits>
#include <map>
#include <vector>

#include <cmath>

#include <catch2/catch.hpp>

#include "chemistry.hpp"
#include "config.hpp"
#include "test_config.hpp"

static void Equilibrate(Cantera::ThermoPhase& thermo,
                        double const temperature,
                        double pressure,
                        Cantera::compositionMap const& composition)
{
  try
    {
      // set thermo state
      thermo.setState_TPX(temperature, pressure, composition);

      // compute equilibrium solution
      thermo.equilibrate("HP");
    }
  catch (Cantera::CanteraError const& err)
    {
      std::cerr << err.what() << std::endl;
    }
}

struct ThermoOutput
{
  double temperature;
  double pressure;
  double density;
  double molar_enthalpy;
  double molar_entropy;
  double molar_specific_heat;
};

static ThermoOutput GetThermoOutput(Cantera::ThermoPhase const& thermo, bool verbose = false)
{
  ThermoOutput output;
  output.temperature = thermo.temperature();
  output.pressure = thermo.pressure();
  output.density = thermo.density();
  output.molar_enthalpy = thermo.enthalpy_mole();
  output.molar_entropy = thermo.entropy_mole();
  output.molar_specific_heat = thermo.cp_mole();

  if (verbose)
    {
      Cantera::writelog("\nThermo properties at the equilibrium state:\n"
                        "Temperature:          {:14.5g} K\n"
                        "Pressure:             {:14.5g} Pa\n"
                        "Density:              {:14.5g} kg/m3\n"
                        "Molar Enthalpy:       {:14.5g} J/kmol\n"
                        "Molar Entropy:        {:14.5g} J/kmol-K\n"
                        "Molar cp:             {:14.5g} J/kmol-K\n",
                        output.temperature,
                        output.pressure,
                        output.density,
                        output.molar_enthalpy,
                        output.molar_entropy,
                        output.molar_specific_heat);
    }

  return output;
}

struct TransportOutput
{
  double viscosity;
  double thermal_conductivity;
};

static TransportOutput GetTransOutput(Cantera::Transport& trans,
                                      bool verbose = false)
{
  TransportOutput output;
  output.viscosity = trans.viscosity();
  output.thermal_conductivity = trans.thermalConductivity();

  if (verbose)
    {
      Cantera::writelog("\nTransport properties at the equilibrium state:\n"
                        "Viscosity:            {:14.5g} kg/m-s\n"
                        "Thermal Conductivity: {:14.5g} W/m-K\n",
                        output.viscosity,
                        output.thermal_conductivity);
    }

  return output;
}

struct RateOfProgress
{
  std::vector<double> fwd;
  std::vector<double> rev;
  std::vector<double> net;
  RateOfProgress()
      : fwd(), rev(), net()
  {
  }
  RateOfProgress(size_t n)
      : fwd(n), rev(n), net(n)
  {
  }
};

struct KineticsOutput
{
  RateOfProgress rop;
};

static KineticsOutput GetKineticsOutput(Cantera::Kinetics& kinetics,
                                        bool verbose = false)
{
  KineticsOutput output;
  size_t const n_rxns = kinetics.nReactions();

  output.rop = RateOfProgress(n_rxns);
  kinetics.getFwdRatesOfProgress(output.rop.fwd.data());
  kinetics.getRevRatesOfProgress(output.rop.rev.data());
  kinetics.getNetRatesOfProgress(output.rop.net.data());

  if (verbose)
    {
      Cantera::writelog("\nReactions and their forward, reverse and net rates of progress:\n");
      for (size_t i = 0; i < n_rxns; i++)
        {
          Cantera::writelog("{:6s} {:35s} {:14.5g} {:14.5g} {:14.5g}  kmol/m3/s\n",
                            'R' + std::to_string(i + 1),
                            kinetics.reactionString(i),
                            output.rop.fwd[i],
                            output.rop.rev[i],
                            output.rop.net[i]);
        }
    }

  return output;
}

static void TestThermoOutput(ThermoOutput const& output)
{
  CHECK(Approx(1.706732e+03).epsilon(TestConfig::tolerance) == output.temperature);
  CHECK(Approx(2.026500e+05).epsilon(TestConfig::tolerance) == output.pressure);
  CHECK(Approx(3.239982e-01).epsilon(TestConfig::tolerance) == output.density);
  CHECK(Approx(-5.626613e+06).epsilon(TestConfig::tolerance) == output.molar_enthalpy);
  CHECK(Approx(2.433695e+05).epsilon(TestConfig::tolerance) == output.molar_entropy);
  CHECK(Approx(3.734765e+04).epsilon(TestConfig::tolerance) == output.molar_specific_heat);
}

static void TestTransportOutput(TransportOutput const& output)
{
  CHECK(Approx(5.837424e-05).epsilon(TestConfig::tolerance) == output.viscosity);
  CHECK(Approx(1.747307e-01).epsilon(TestConfig::tolerance) == output.thermal_conductivity);
}

static void TestKineticsOutput(KineticsOutput const& output)
{
  RateOfProgress const& rop = output.rop;
  for (size_t i = 0; i < rop.net.size(); i++)
    {
      CHECK_THAT(output.rop.net[i], Catch::Matchers::WithinAbs(0.0, TestConfig::tolerance));
    }
}

SCENARIO("Chemical equilibtium state can be computed", "[chemistry]")
{

  GIVEN("Chemistry parameters")
  {
    std::string const mechanism = "../data/gri30.cti";
    double const temperature = 500.0;
    double const pressure = 2.0 * Cantera::OneAtm;
    Cantera::compositionMap const composition = {
        {"CH4", 1.00},
        {"O2", 1.00},
        {"N2", 3.76}};

    try
      {

        WHEN("Thermo is requested")
        {

          Chemistry::Thermo chemistry(mechanism);
          Cantera::ThermoPhase& thermo = chemistry.thermo();
          Equilibrate(thermo, temperature, pressure, composition);
          THEN("Check thermo properties")
          {
            REQUIRE_NOTHROW(TestThermoOutput(GetThermoOutput(thermo, TestConfig::verbose)));
          }
        }

        WHEN("Thermoy and transport are requested")
        {
          Chemistry::TransportModel const transport_model = Chemistry::TransportModel::MixtureAveraged;
          Chemistry::ThermoTransport chemistry(mechanism, transport_model);
          Cantera::ThermoPhase& thermo = chemistry.thermo();
          Equilibrate(thermo, temperature, pressure, composition);

          // test output
          THEN("Check thermo and transport properties")
          {
            REQUIRE_NOTHROW(TestThermoOutput(GetThermoOutput(thermo, TestConfig::verbose)));
            REQUIRE_NOTHROW(TestTransportOutput(GetTransOutput(chemistry.transport(), TestConfig::verbose)));
          }
        }

        WHEN("Thermo and kinetics and kinetics are requested")
        {
          Chemistry::ThermoKinetics chemistry(mechanism);
          Cantera::ThermoPhase& thermo = chemistry.thermo();
          Equilibrate(thermo, temperature, pressure, composition);

          // test output
          THEN("Check thermo and kinetic properties")
          {
            REQUIRE_NOTHROW(TestThermoOutput(GetThermoOutput(thermo, TestConfig::verbose)));
            REQUIRE_NOTHROW(TestKineticsOutput(GetKineticsOutput(chemistry.kinetics(), TestConfig::verbose)));
          }
        }

        WHEN("Thermo, transport and kinetics are requested")
        {
          Chemistry::TransportModel const transport_model = Chemistry::TransportModel::MixtureAveraged;
          Chemistry::ThermoTransportKinetics chemistry(mechanism, transport_model);
          Cantera::ThermoPhase& thermo = chemistry.thermo();
          Equilibrate(thermo, temperature, pressure, composition);

          // test output
          THEN("Check thermo, transport and kinetic properties")
          {
            REQUIRE_NOTHROW(TestThermoOutput(GetThermoOutput(thermo, TestConfig::verbose)));
            REQUIRE_NOTHROW(TestTransportOutput(GetTransOutput(chemistry.transport(), TestConfig::verbose)));
            REQUIRE_NOTHROW(TestKineticsOutput(GetKineticsOutput(chemistry.kinetics(), TestConfig::verbose)));
          }
        }
      }
    catch (Cantera::CanteraError const& err)
      {
        std::cout << err.what() << std::endl;
      }
  }
}

TEST_CASE("Chemistry objects can be moved", "[chemistry]")
{

  enum
  {
    GRI12 = 0,
    GRI211 = 1,
    GRI30 = 2
  };
  std::map<int, std::string> const mechanism = {
      std::make_pair(GRI12, "../data/gri12.cti"),
      std::make_pair(GRI211, "../data/gri211.cti"),
      std::make_pair(GRI30, "../data/gri30.cti"),
  };
  double const temperature = 500.0;
  double const pressure = 2.0 * Cantera::OneAtm;
  Cantera::compositionMap const composition = {
      {"CH4", 1.00},
      {"O2", 1.00},
      {"N2", 3.76}};

  std::vector<Chemistry::ThermoTransportKinetics> chemistry;
  chemistry.reserve(mechanism.size());

  Chemistry::TransportModel const transport_model = Chemistry::TransportModel::MixtureAveraged;

  // emplace back first instance GRI-Mec 1.2
  chemistry.emplace_back(mechanism.at(GRI12), transport_model);
  {
    {
      Cantera::ThermoPhase const& thermo = chemistry[GRI12].thermo();
      CHECK(32 == thermo.nSpecies());
      CHECK(13 == thermo.speciesIndex("CH4"));
    }
    {
      Cantera::Kinetics const& kinetics = chemistry[GRI12].kinetics();
      CHECK(177 == kinetics.nReactions());
      CHECK("C2H4 + CH3 <=> C2H3 + CH4" == kinetics.reactionString(163));
    }
  }

  // push back second instance GRI-Mec 2.11
  chemistry.push_back(std::move(Chemistry::ThermoTransportKinetics(mechanism.at(GRI211), transport_model)));
  {
      {Cantera::ThermoPhase const& thermo = chemistry[GRI211].thermo();
  CHECK(49 == thermo.nSpecies());
  CHECK(6 == thermo.speciesIndex("HO2"));
}
{
  Cantera::Kinetics const& kinetics = chemistry[GRI211].kinetics();
  CHECK(279 == kinetics.nReactions());
  CHECK("CH3O + H <=> CH3 + OH" == kinetics.reactionString(65));
}
}

// insert third instance at the end of the vecor GRI-Mec 3.0
{
  Chemistry::ThermoTransportKinetics instance_1;
  instance_1.Create(mechanism.at(GRI30), transport_model);
  Chemistry::ThermoTransportKinetics instance_2;
  instance_2 = std::move(instance_1);
  Chemistry::ThermoTransportKinetics instance_3(std::move(instance_2));
  chemistry.insert(chemistry.end(), std::move(instance_3));
}

{
  Chemistry::ThermoTransportKinetics& chem_gri30 = chemistry.back();
  Equilibrate(chem_gri30.thermo(), temperature, pressure, composition);
  REQUIRE_NOTHROW(TestThermoOutput(GetThermoOutput(chem_gri30.thermo())));
  REQUIRE_NOTHROW(TestTransportOutput(GetTransOutput(chem_gri30.transport())));
  REQUIRE_NOTHROW(TestKineticsOutput(GetKineticsOutput(chem_gri30.kinetics())));
}
}
