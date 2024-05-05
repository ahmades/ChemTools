// Adapted from:
// https://www.codeproject.com/Articles/820116/Embedding-Python-program-in-a-C-Cplusplus-code

#ifndef EMBEDDED_PYTHON_HELPER_HPP
#define EMBEDDED_PYTHON_HELPER_HPP

#include <Python.h>

namespace CppPy
{

  /*
     Initialises and finalises the Python interprete
     When needed, a single instance should be created during initialisation
  */

  class Instance
  {

  public:
    Instance()
    {
      Py_Initialize();
    }

    ~Instance()
    {
      Py_Finalize();
    }
  };

  // Helper class for PyObject management

  class Object
  {

  private:
    PyObject* py_obj;

  public:
    Object()
        : py_obj(nullptr)
    {
    }

    Object(PyObject* const py_obj_)
        : py_obj(py_obj_)
    {
    }

    ~Object()
    {
      Release();
    }

    operator PyObject*()
    {
      return py_obj;
    }

    PyObject* operator=(PyObject* const py_obj_)
    {
      py_obj = py_obj_;
      return py_obj;
    }

    void IncreaseReference()
    {
      if (py_obj)
        {
          Py_INCREF(py_obj);
        }
    }

    void DecreaseReference()
    {
      if (py_obj)
        {
          Py_DECREF(py_obj);
        }
    }

    void Release()
    {
      DecreaseReference();
      py_obj = nullptr;
    }

    bool IsValid()
    {
      return py_obj ? true : false;
    }
  };

} // namespace CppPy

#endif // EMBEDDED_PYTHON_HELPER_HPP
