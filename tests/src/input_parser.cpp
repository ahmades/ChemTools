#include <algorithm>
#include <string>
#include <vector>

#include <catch2/catch.hpp>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/printf.h>

#include <boost/filesystem.hpp>

#include "yaml-cpp/yaml.h"

#include "units.hpp"

#include <catch2/catch.hpp>

#include "input/parser.hpp"
#include "input/types.hpp"
#include "test_config.hpp"

TEST_CASE("Input can be parsed", "[input]")
{

  // inexistent config
  {
    std::string const config = "../non_existent_config.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), YAML::BadFile);
  }

  // invalid node: pressure value missing
  {
    std::string const config = "../tests/data/yaml_input/config_invalid_node.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), YAML::InvalidNode);
  }

  // invalid conversion: pressure value is a string instead of double
  {
    std::string const config = "../tests/data/yaml_input/config_bad_conversion.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), YAML::BadConversion);
  }

  // invalid mechanism: wrong path
  {
    std::string const config = "../tests/data/yaml_input/config_invalid_mech.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), boost::filesystem::filesystem_error);
  }

  // invalid value: negative pressure
  {
    std::string const config = "../tests/data/yaml_input/config_negative_press.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // invalid value: negative temperature
  {
    std::string const config = "../tests/data/yaml_input/config_negative_temp.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // duplicate species in composition tuples
  {
    std::string const config = "../tests/data/yaml_input/config_duplicate_spec.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // negative species in composition tuples
  {
    std::string const config = "../tests/data/yaml_input/config_negative_spec.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // invalid species in composition tuples
  {
    std::string const config = "../tests/data/yaml_input/config_invalid_spec.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // nonconservative composition
  {
    std::string const config = "../tests/data/yaml_input/config_nonconservative_comp.yaml";
    Input::Parser parser;
    REQUIRE_THROWS_AS(parser = Input::Parser(config), std::runtime_error);
  }

  // different types of case sets with variable number of cases
  {
    std::string const config = "../tests/data/yaml_input/config_multiple_reactors.yaml";
    Input::Parser parser;
    REQUIRE_NOTHROW(parser = Input::Parser(config));
    Input::Parameters const& parameters = parser.GetParameters();
    Input::Reactor const& reactor = parameters.reactor;

    // const pressure
    {
      std::vector<Input::CaseSet> const& case_sets = reactor.const_P;
      REQUIRE(2 == case_sets.size());
      {
        std::vector<Input::Case> const& cases = case_sets[0].Cases();
        REQUIRE(3 == cases.size());
      }
      {
        std::vector<Input::Case> const& cases = case_sets[1].Cases();
        REQUIRE(2 == cases.size());
      }
    }

    // const volunme
    {
      std::vector<Input::CaseSet> const& case_sets = reactor.const_V;
      REQUIRE(1 == case_sets.size());
      std::vector<Input::Case> const& cases = case_sets[0].Cases();
      REQUIRE(2 == cases.size());
    }

    // const temperature and pressure
    {
      std::vector<Input::CaseSet> const& case_sets = reactor.const_TP;
      REQUIRE(1 == case_sets.size());
      std::vector<Input::Case> const& cases = case_sets[0].Cases();
      REQUIRE(1 == cases.size());
    }

    // const temperature and volume
    {
      std::vector<Input::CaseSet> const& case_sets = reactor.const_TV;
      REQUIRE(1 == case_sets.size());
      std::vector<Input::Case> const& cases = case_sets[0].Cases();
      REQUIRE(2 == cases.size());
    }
  }

  {
    std::string const config = "../tests/data/yaml_input/config.yaml";
    Input::Parser parser;
    REQUIRE_NOTHROW(parser = Input::Parser(config));

    Input::Parameters const& parameters = parser.GetParameters();
    Input::Reactor const& reactor = parameters.reactor;

    REQUIRE_FALSE(reactor.const_P.empty());
    REQUIRE(reactor.const_V.empty());
    REQUIRE(reactor.const_TP.empty());
    REQUIRE(reactor.const_TV.empty());

    // const pressure reactor tests

    std::vector<Input::CaseSet> const& case_sets = reactor.const_P;
    REQUIRE(1 == case_sets.size());

    std::vector<Input::Case> const& cases = case_sets[0].Cases();
    REQUIRE(4 == cases.size());

    // case 0
    {
      Input::Case const& cas = cases[0];
      {
        Input::RealScalar const& pressure = cas.pressure;
        REQUIRE(101325.0 == pressure.Value());
        REQUIRE(units::precise::pascal == pressure.InputUnit());
        REQUIRE(units::precise::pascal == pressure.Unit());
      }

      {
        Input::RealScalar const& temperature = cas.temperature;
        REQUIRE(1300.0 == temperature.Value());
        REQUIRE(units::precise::Kelvin == temperature.InputUnit());
        REQUIRE(units::precise::Kelvin == temperature.Unit());
      }

      {
        Input::Composition const& composition = cas.composition;
        REQUIRE(Input::CompositionInputType::MoleFraction == composition.Type());
        Input::CompositionPairs pairs = composition.Pairs();
        REQUIRE(3 == pairs.size());
        double sum = 0.0;
        for (size_t i = 0; i < pairs.size(); ++i)
          {
            if (0 == i)
              {
                REQUIRE("CH4" == pairs[i].first);
                REQUIRE(0.1 == pairs[i].second);
              }
            if (1 == i)
              {
                REQUIRE("O2" == pairs[i].first);
                REQUIRE(0.5 == pairs[i].second);
              }
            if (2 == i)
              {
                REQUIRE("N2" == pairs[i].first);
                REQUIRE(0.4 == pairs[i].second);
              }
            sum += pairs[i].second;
          }
        REQUIRE(1.0 == sum);
      }

      {
        Input::RealScalar const& time = cas.time;
        REQUIRE(1.0e-3 == time.Value());
        REQUIRE(units::precise::ms == time.InputUnit());
        REQUIRE(units::precise::s == time.Unit());
      }
    }

    // case 1
    {
      Input::Case const& cas = cases[1];
      {
        Input::RealScalar const& pressure = cas.pressure;
        REQUIRE(202650.0 == pressure.Value());
        REQUIRE(units::precise::pressure::atm == pressure.InputUnit());
        REQUIRE(units::precise::pascal == pressure.Unit());
      }

      {
        Input::RealScalar const& temperature = cas.temperature;
        REQUIRE(1298.15 == temperature.Value());
        REQUIRE(units::precise::degC == temperature.InputUnit());
        REQUIRE(units::precise::Kelvin == temperature.Unit());
      }

      {
        Input::Composition const& composition = cas.composition;
        REQUIRE(Input::CompositionInputType::Concentration == composition.Type());
        REQUIRE(units::precise_unit(units::milli * units::mol / units::L) == composition.InputUnit());
        REQUIRE(units::precise_unit(units::mol / units::m.pow(3)) == composition.Unit());
        Input::CompositionPairs pairs = composition.Pairs();
        REQUIRE(4 == pairs.size());
        for (size_t i = 0; i < pairs.size(); ++i)
          {
            if (0 == i)
              {
                REQUIRE("CH4" == pairs[i].first);
                REQUIRE(0.1 == pairs[i].second);
              }
            if (1 == i)
              {
                REQUIRE("O2" == pairs[i].first);
                REQUIRE(0.5 == pairs[i].second);
              }
            if (2 == i)
              {
                REQUIRE("N2" == pairs[i].first);
                REQUIRE(0.35 == pairs[i].second);
              }
            if (3 == i)
              {
                REQUIRE("H2" == pairs[i].first);
                REQUIRE(0.05 == pairs[i].second);
              }
          }
      }

      {
        Input::RealScalar const& time = cas.time;
        REQUIRE(60.0 == time.Value());
        REQUIRE(units::precise::min == time.InputUnit());
        REQUIRE(units::precise::s == time.Unit());
      }
    }

    // case 2
    {
      Input::Case const& cas = cases[2];
      {
        Input::RealScalar const& pressure = cas.pressure;
        REQUIRE(150000.0 == pressure.Value());
        REQUIRE(units::precise::bar == pressure.InputUnit());
        REQUIRE(units::precise::pascal == pressure.Unit());
      }

      {
        Input::RealScalar const& temperature = cas.temperature;
        REQUIRE(1368.15 == temperature.Value());
        REQUIRE(units::precise::degF == temperature.InputUnit());
        REQUIRE(units::precise::Kelvin == temperature.Unit());
      }

      {
        Input::Composition const& composition = cas.composition;
        REQUIRE(Input::CompositionInputType::Concentration == composition.Type());
        REQUIRE(units::precise_unit(units::milli * units::mol / units::L) == composition.InputUnit());
        REQUIRE(units::precise_unit(units::mol / units::m.pow(3)) == composition.Unit());
        Input::CompositionPairs pairs = composition.Pairs();
        REQUIRE(4 == pairs.size());
        for (size_t i = 0; i < pairs.size(); ++i)
          {
            if (0 == i)
              {
                REQUIRE("CH3" == pairs[i].first);
                REQUIRE(0.03 == pairs[i].second);
              }
            if (1 == i)
              {
                REQUIRE("CH4" == pairs[i].first);
                REQUIRE(0.07 == pairs[i].second);
              }
            if (2 == i)
              {
                REQUIRE("O2" == pairs[i].first);
                REQUIRE(0.5 == pairs[i].second);
              }
            if (3 == i)
              {
                REQUIRE("N2" == pairs[i].first);
                REQUIRE(0.4 == pairs[i].second);
              }
          }
      }

      {
        Input::RealScalar const& time = cas.time;
        REQUIRE(3600.0 == time.Value());
        REQUIRE(units::precise::time::h == time.InputUnit());
        REQUIRE(units::precise::s == time.Unit());
      }
    }

    // case 3
    {
      Input::Case const& cas = cases[3];
      {
        Input::RealScalar const& pressure = cas.pressure;
        REQUIRE(201326.0 == pressure.Value());
        REQUIRE(units::precise::pascal == pressure.InputUnit());
        REQUIRE(units::precise::pascal == pressure.Unit());
      }

      {
        Input::RealScalar const& temperature = cas.temperature;
        REQUIRE(1250.0 == temperature.Value());
        REQUIRE(units::precise::Kelvin == temperature.InputUnit());
        REQUIRE(units::precise::Kelvin == temperature.Unit());
      }

      {
        Input::Composition const& composition = cas.composition;
        REQUIRE(Input::CompositionInputType::MassFraction == composition.Type());
        Input::CompositionPairs pairs = composition.Pairs();
        REQUIRE(3 == pairs.size());
        double sum = 0.0;
        for (size_t i = 0; i < pairs.size(); ++i)
          {
            if (0 == i)
              {
                REQUIRE("CH4" == pairs[i].first);
                REQUIRE(0.15 == pairs[i].second);
              }
            if (1 == i)
              {
                REQUIRE("O2" == pairs[i].first);
                REQUIRE(0.45 == pairs[i].second);
              }
            if (2 == i)
              {
                REQUIRE("N2" == pairs[i].first);
                REQUIRE(0.4 == pairs[i].second);
              }
            sum += pairs[i].second;
          }
        REQUIRE(1.0 == sum);
      }

      {
        Input::RealScalar const& time = cas.time;
        REQUIRE(1.0 == time.Value());
        REQUIRE(units::precise::s == time.InputUnit());
        REQUIRE(units::precise::s == time.Unit());
      }
    }

    if (TestConfig::verbose)
      {

        for (size_t j = 0; j < case_sets.size(); ++j)
          {
            Input::CaseSet const& case_set = case_sets[j];

            fmt::print(fmt::emphasis::underline | fmt::emphasis::bold,
                       "Set {}\n\n",
                       j + 1);

            fmt::print("Mechanism: Path        = {}\n"
                       "           Description = {}\n\n",
                       case_set.GetMechanism().Path(),
                       case_set.GetMechanism().Description());

            if (case_set.DescriptionIsDefined())
              {
                fmt::print("Description: {}\n\n", case_set.Description());
              }

            for (auto cas : case_set.Cases())
              {
                std::string species_str;
                std::for_each(cas.composition.Pairs().begin(),
                              cas.composition.Pairs().end(),
                              [&species_str](const auto& pair) { species_str += fmt::format("{}:{} ", pair.first, pair.second); });
                if (cas.DescriptionIsDefined())
                  {
                    fmt::print("Case:\n      Description: {}\n", cas.Description());
                  }
                else
                  {
                    fmt::print("Case:\n");
                  }

                fmt::print("      P = {} {}\n"
                           "      T = {} {}\n"
                           "      Y = {} {}\n"
                           "      t = {} {}\n\n",
                           cas.pressure.Value(),
                           units::to_string(cas.pressure.Unit()),
                           cas.temperature.Value(),
                           units::to_string(cas.temperature.Unit()),
                           species_str,
                           units::to_string(cas.composition.Unit()),
                           cas.time.Value(),
                           units::to_string(cas.time.Unit()));
              }
          }
      }
  }
}
