#include <fmt/format.h>
#include <fmt/printf.h>

#include <catch2/catch.hpp>

#include <boost/filesystem.hpp>

#include "embedded_python/functions.hpp"
#include "test_config.hpp"

TEST_CASE("Python function with one parameter and returning a pair can be called from c++", "[embedde_python]")
{

  std::string const module_path("../tests/data/embedded_python/script.py");
  std::string const function_name("function_1");
  CppPy::FunctionP1R2 function;
  REQUIRE_NOTHROW(function.Load(module_path, function_name));

  double x{2.0};
  double x2{0.0};
  double x3{0.0};
  REQUIRE(true == function.Evaluate(x, x2, x3));
  REQUIRE(x * x == x2);
  REQUIRE(x * x2 == x3);

  if (TestConfig::verbose)
    {
      fmt::print("\nIf x = {} then x^2 = {} and x^3 = {}\n", x, x2, x3);
    }
}
