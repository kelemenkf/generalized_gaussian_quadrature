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
        py::function,
        std::function<double(const double&)>,
        std::function<double(const double&, const double&)>,
        std::function<double(const double&, const double&, const double&)>,
        std::function<double(const double&, const double&, const double&, const double&)>
    >;


    InputFunction functionVariant;
    std::vector<std::vector<double>> paramSpace; 
    size_t numberOfParameters;


public:
    bool pythonFlag;


public:
    // FunctionHandler(py::function inputFunction, Parameter... parameters) : function(inputFunction)
    // {
    //     (paramSpace.push_back(parameters),...);
    //     numberOfParameters = sizeof...(parameters);
    //     pythonFlag = true;
    // }
    
    
    FunctionHandler(InputFunction inputFunction, Parameter... parameters) : functionVariant(inputFunction) 
    {
        (paramSpace.push_back(parameters),...);
        numberOfParameters = sizeof...(parameters);
        pythonFlag = false;
        if (functionVariant.index() == 0)
            pythonFlag = true;
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    double callFunction(double x, double param1 = 10, double param2 = 0.5, double param3 = 0) 
    {
        size_t index = functionVariant.index(); 
        if (index == 1 && numberOfParameters == 0) 
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&)>>(&functionVariant)) 
            {
                return (*functionPtr)(x);
            }
        } 
        else if (index == 2 && numberOfParameters == 1)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1);
            }
        }
        else if (index == 3 && numberOfParameters == 2)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2);
            }
        }
        else if (index == 4 && numberOfParameters == 3)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double& , const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2, param3);
            }
        }
        else if (index == 0)
        {
            callFunctionPython(x, param1, param2, param3);
        }
        else
        {
            throw std::runtime_error("Number of parameters is not consistent with passed function");
        }
    }


    double callFunctionPython(double x, double param1 = 10, double param2 = 0.5, double param3 = 0) 
    {
        auto function = std::get_if<py::function>(&functionVariant);
        if (numberOfParameters == 0) 
        {
            return (*function)(x).template cast<double>();  
        } 
        else if (numberOfParameters == 1)
        {
            return (*function)(x, param1).template cast<double>();
        }
        else if (numberOfParameters == 2)
        {
            return (*function)(x, param1, param2).template cast<double>();
        }
        else if (numberOfParameters == 3)
        {
            return (*function)(x, param1, param2, param3).template cast<double>();
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