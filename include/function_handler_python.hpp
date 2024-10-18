#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdarg> 
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

template<typename... Parameter>
class FunctionHandlerPython
{
private:
    std::vector<std::vector<double>> paramSpace; 
    size_t numberOfParameters;
    py::function function;


public:
    FunctionHandler(py::function inputFunction, Parameter... parameters) : function(inputFunction)
    {
        (paramSpace.push_back(parameters),...);
        numberOfParameters = sizeof...(parameters);
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    double callFunctionPython(double x, double param1 = 10, double param2 = 0.5, double param3 = 0) 
    {
        if (numberOfParameters == 0) 
        {
            return function(x).template cast<double>();  
        } 
        else if (numberOfParameters == 1)
        {
            return function(x, param1).template cast<double>();
        }
        else if (numberOfParameters == 2)
        {
            return function(x, param1, param2).template cast<double>();
        }
        else if (numberOfParameters == 3)
        {
            return function(x, param1, param2, param3).template cast<double>();
        }
        else
        {
            throw std::runtime_error("Number of parameters is not consistent with passed function");
        }
    }
    

    double getNumberOfParameters()
    {
        return numberOfParameters;
    }
};