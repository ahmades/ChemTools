#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/system/error_code.hpp>

#include "embedded_python/functions.hpp"

namespace CppPy
{

  FunctionP1R2::FunctionP1R2()
      : function(nullptr), arg(nullptr), value(nullptr)
  {
  }

  FunctionP1R2::FunctionP1R2(std::string const& module_path, std::string const& function_name)
  {
    Load(module_path, function_name);
  }

  void FunctionP1R2::Load(std::string const& module_path, std::string const& function_name)
  {

    if (!Py_IsInitialized())
      {
        throw std::runtime_error("FunctionP1R2 error: "
                                 "Python interpreter has not been initialised.");
      }

    // check if module exists
    boost::filesystem::path const bfsp(module_path);
    boost::system::error_code error_code;
    if (!boost::filesystem::exists(bfsp, error_code))
      {
        BOOST_FILESYSTEM_THROW(boost::filesystem::filesystem_error("FunctionP1R2 error",
                                                                   bfsp,
                                                                   error_code));
      }

    // add module parent path to system path
    // get system path
    CppPy::Object sys_path = PySys_GetObject("path");
    if (!sys_path.IsValid() || !PyList_Check(sys_path))
      {
        throw std::runtime_error("FunctionP1R2 error: "
                                 "Could not access python sys.path.");
      }

    // get module parent path
    CppPy::Object module_parent_path = PyUnicode_FromString(bfsp.parent_path().string().c_str());

    // append module parent path to system path
    if (PyList_Append(sys_path, module_parent_path) != 0)
      {
        throw std::runtime_error("FunctionP1R2 error: "
                                 "Could not append module path to system path.");
      }

    // import module
    CppPy::Object module = PyImport_ImportModule(bfsp.filename().stem().string().c_str());

    // get function from module
    if (module.IsValid())
      {
        function = PyObject_GetAttrString(module, function_name.c_str());
        if (!function.IsValid() || !PyCallable_Check(function))
          {
            throw std::invalid_argument("FunctionP1R2 error: Could not find function \"" +
                                        function_name +
                                        "\" in imported module \"" +
                                        module_path +
                                        "\".");
          }
      }
    else
      {
        throw std::invalid_argument("FunctionP1R2 error: Could not import module \"" +
                                    module_path +
                                    "\".");
      }
  }

  bool FunctionP1R2::Evaluate(double const param,
                              double& val_1,
                              double& val_2)
  {
    arg = Py_BuildValue("(d)", param);
    if (!arg.IsValid())
      {
        return false;
      }
    value = PyObject_CallObject(function, arg);
    if (!value.IsValid())
      {
        return false;
      }
    return PyArg_ParseTuple(value, "d|d", &val_1, &val_2);
  }

} // namespace CppPy
