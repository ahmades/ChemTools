#ifndef INPUT_PARSER_HPP
#define INPUT_PARSER_HPP

#include "input/types.hpp"

namespace Input
{

  class Parser
  {
  private:
    std::string config_file;
    Parameters parameters;

    // YAML parser
    void Parse();

    // vaildation of parsed input
    void Validate() const;

  public:
    // constructors
    Parser();
    Parser(std::string const& config_file_);

    // parameters getter
    Parameters const& GetParameters() const;
  };

} // namespace Input

#endif // INPUT_PARSER_HPP
