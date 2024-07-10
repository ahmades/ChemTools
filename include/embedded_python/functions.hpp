#include <Python.h>

#include "embedded_python/helper.hpp"

namespace CppPy
{

  /*
     FunctionP1R2 loads and evaluates a python function
     taking 1 parameter (P1) and returning 2 values (R2).
     The python function msut be defined as:
     def function( param ):
         assign( val_1 )
         assign( val_2 )
         return ( val_1, val_2 )
  */
  class FunctionP1R2
  {

  private:
    Object function;
    Object arg;
    Object value;

  public:
    FunctionP1R2();

    FunctionP1R2(std::string const& module_path,
                 std::string const& function_name);

    void Load(std::string const& module_path,
              std::string const& function_name);

    bool Evaluate(double const param,
                  double& val_1,
                  double& val_2);
  };

} // namespace CppPy
