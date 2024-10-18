#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdarg> 
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

namespace py = pybind11;

template<typename... Parameter>
class FunctionHandler
{
private:
    using InputFunction = std::variant
    <
        std::function<double(const double&)>,
        std::function<double(const double&, const double&)>,
        std::function<double(const double&, const double&, const double&)>,
        std::function<double(const double&, const double&, const double&, const double&)>
    >;


    InputFunction functionVariant;
    std::vector<std::vector<double>> paramSpace; 
    size_t numberOfParameters;
    size_t index; 
    py::function function;


public:
    bool pythonFlag;


public:
    FunctionHandler(py::function inputFunction, Parameter... parameters) : function(inputFunction)
    {
        (paramSpace.push_back(parameters),...);
        numberOfParameters = sizeof...(parameters);
        pythonFlag = true;
    }
    
    
    FunctionHandler(InputFunction inputFunction, Parameter... parameters) : functionVariant(inputFunction) 
    {
        (paramSpace.push_back(parameters),...);
        numberOfParameters = sizeof...(parameters);
        index = functionVariant.index();
        pythonFlag = false;
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    double callFunction(double x, double param1 = 10, double param2 = 0.5, double param3 = 0) 
    {
        if (index == 0 && numberOfParameters == 0) 
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&)>>(&functionVariant)) 
            {
                return (*functionPtr)(x);
            }
        } 
        else if (index == 1 && numberOfParameters == 1)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1);
            }
        }
        else if (index == 2 && numberOfParameters == 2)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2);
            }
        }
        else if (index == 3 && numberOfParameters == 3)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double& , const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2, param3);
            }
        }
        else
        {
            throw std::runtime_error("Number of parameters is not consistent with passed function");
        }
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


    double evaluateFunctionOnParameterSpace()
    {

    }


    double getNumberOfParameters()
    {
        return numberOfParameters;
    }
};