#ifndef INPUT_TYPES_HPP
#define INPUT_TYPES_HPP

#include <boost/bimap.hpp>
#include <functional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>

#include "units.hpp"

namespace Input
{

  class DescriptionBase
  {
  private:
    std::string description;
    bool is_defined;

  public:
    DescriptionBase()
        : description(), is_defined{false}
    {
    }

    // setters
    void SetDescription(std::string const& description_)
    {
      description = description_;
      is_defined = true;
    }

    // getters
    std::string const& Description() const
    {
      return description;
    }

    bool DescriptionIsDefined() const
    {
      return is_defined;
    }
  };

  class Mechanism : public DescriptionBase
  {
  private:
    std::string path;

  public:
    // constructor: can be omitted, default ctor is OK
    Mechanism()
        : path()
    {
    }

    // setters
    void SetPath(std::string const& path_)
    {
      path = path_;
    }

    // getters
    std::string const& Path() const
    {
      return path;
    }
  };

  class RealScalar
  {
  private:
    units::precise_unit inp_unit;
    units::precise_unit unit;
    std::string unit_str;
    double value;

  public:
    // constructor
    RealScalar(units::precise_unit const unit_)
        : inp_unit{units::precise::invalid}, unit{unit_}, value{std::numeric_limits<double>::signaling_NaN()}
    {
      if (!units::is_valid(unit))
        {
          throw(std::invalid_argument("Invalid unit."));
        }
      unit_str = units::to_string(unit);
    }

    // setters
    void Set(double const inp_value, std::string const& inp_unit_str)
    {
      inp_unit = units::unit_from_string(inp_unit_str);
      if (!units::is_valid(inp_unit))
        {
          throw(std::invalid_argument(inp_unit_str
                                      + " is not a valid unit."));
        }
      if (!inp_unit.is_convertible(unit))
        {
          throw(std::invalid_argument(inp_unit_str
                                      + " is not convertible to "
                                      + unit_str + '.'));
        }
      value = units::convert(inp_value, inp_unit, unit);
    }

    // getters
    double Value() const
    {
      return value;
    }

    units::precise_unit const Unit() const
    {
      return unit;
    }

    units::precise_unit const InputUnit() const
    {
      return inp_unit;
    }
  };

  using CompositionPairs = std::vector<std::pair<std::string, double>>;

  enum class CompositionInputType
  {
    Invalid,
    Concentration,
    MassFraction,
    MoleFraction
  };

  class Composition
  {
  private:
    CompositionInputType type;
    units::precise_unit inp_unit;
    units::precise_unit unit;
    // maybe change this to 2 vectors: names and values?
    CompositionPairs composition_pairs;

  public:
    Composition()
        : type{CompositionInputType::Invalid}, inp_unit{units::precise::invalid}, unit{units::precise::invalid}, composition_pairs()
    {
    }

    void Set(CompositionPairs const& composition_pairs_, std::string const& type_str_, std::string const& inp_unit_str_ = "")
    {
      composition_pairs = composition_pairs_;
      std::string const type_str = boost::algorithm::to_lower_copy(type_str_);
      if ("concentration" == type_str)
        {
          // only dimensional case
          type = CompositionInputType::Concentration;

          inp_unit = units::unit_from_string(inp_unit_str_);
          if (!units::is_valid(inp_unit))
            {
              throw(std::invalid_argument(inp_unit_str_
                                          + " is not a valid unit"));
            }

          unit = units::precise_unit(units::mol / units::m.pow(3));
          if (!inp_unit.is_convertible(unit))
            {
              throw(std::invalid_argument(inp_unit_str_
                                          + " is not convertible to "
                                          + units::to_string(unit)));
            }

          // convert to internal units
          double const conv_factor = units::convert(1.0, inp_unit, unit);
          for (auto& pair : composition_pairs)
            {
              pair.second *= conv_factor;
            }
        }
      else if ("mass fraction" == type_str
               || "mole fraction" == type_str)
        {
          // dimensionless cases: ratio, no conv needed
          type = ("mass fraction" == type_str)
                     ? CompositionInputType::MassFraction
                     : CompositionInputType::MoleFraction;
          unit = inp_unit = units::precise::ratio;
        }
      else
        {
          throw(std::invalid_argument("Composition type is set to \""
                                      + type_str
                                      + "\". The valid types are: \"concentration\", "
                                        "\"mass fraction\" or \"mole fraction\"."));
        }
    }

    CompositionPairs const& Pairs() const
    {
      return composition_pairs;
    }

    std::vector<std::string> const Species() const
    {
      std::vector<std::string> species;
      for (auto const& pair : composition_pairs)
        {
          species.push_back(pair.first);
        }
      return species;
    }

    std::vector<double> const Values() const
    {
      std::vector<double> values;
      for (auto const& pair : composition_pairs)
        {
          values.push_back(pair.second);
        }
      return values;
    }

    units::precise_unit const Unit() const
    {
      return unit;
    }

    units::precise_unit const InputUnit() const
    {
      return inp_unit;
    }

    CompositionInputType Type() const
    {
      return type;
    }
  };

  class Case : public DescriptionBase
  {
  public:
    RealScalar pressure;
    RealScalar temperature;
    Composition composition;
    RealScalar time;

    // constructor
    Case()
        : pressure(units::precise::pascal), temperature(units::precise::Kelvin), composition(), time{units::precise::second}
    {
    }
  };

  class CaseSetBase : public DescriptionBase
  {
  private:
    Mechanism mechanism;

  public:
    // constructor: can be omitted, default ctor is OK
    CaseSetBase()
        : mechanism()
    {
    }

    // setters
    void SetMechanism(Mechanism const& mechanism_)
    {
      mechanism = mechanism_;
    }

    // getters
    Mechanism const& GetMechanism() const
    {
      return mechanism;
    }
  };

  class CaseSet : public CaseSetBase
  {
  private:
    std::vector<Case> cases;

  public:
    // constructor: can be omitted, default ctor is OK
    CaseSet()
        : CaseSetBase(), cases()
    {
    }

    // setters
    void AddCase(Case const& cas)
    {
      cases.push_back(cas);
    }

    // getters
    std::vector<Case> const& Cases() const
    {
      return cases;
    }

    Case const& GetCase(size_t i) const
    {
      return cases[i];
    }

    Case const& operator[](size_t i)
    {
      return cases[i];
    }
  };

  class Reactor
  {
  public:
    std::vector<Input::CaseSet> const_P;  // const pressure
    std::vector<Input::CaseSet> const_V;  // const volume
    std::vector<Input::CaseSet> const_TP; // const temperature and pressure
    std::vector<Input::CaseSet> const_TV; // const temperature and volume

    Reactor()
        : const_P(), const_V(), const_TP(), const_TV()
    {
    }
  };

  class Parameters
  {
  public:
    Reactor reactor;

    Parameters()
        : reactor()
    {
    }
  };

} // namespace Input

#endif // INPUT_TYPES_HPP
