#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdarg> 
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include "utils.hpp"

namespace py = pybind11;

template<typename... Parameter>
class FunctionHandler
{
public: 
    static size_t combinationIndex;


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
    size_t index; 
    std::vector<std::vector<double>> parameterCombinations;
    bool cartesianProduct; 


public:
    FunctionHandler(InputFunction inputFunction, bool cartesianProductInput, Parameter... parameters) 
    : functionVariant(inputFunction), 
    cartesianProduct(cartesianProductInput) 
    {
        if (cartesianProduct)
        {
            index = functionVariant.index();
            (paramSpace.push_back(parameters),...);
            numberOfParameters = sizeof...(parameters);
            buildParameterCombinations();
        }
    }

    
    FunctionHandler(InputFunction inputFunction, bool cartesianProductInput, std::vector<std::vector<double>> paramSpaceInput)
    : paramSpace(paramSpaceInput), cartesianProduct(cartesianProductInput), functionVariant(inputFunction)
    {
        index = functionVariant.index();
        std::cout << "Function index " << index << std::endl;
        buildParameterCombinations();
    }

    
    static void incrementCombinationIndex()
    {
        ++combinationIndex;
    }


    static void resetCombinationIndex()
    {
        combinationIndex = 0;
    }


    static size_t getCombinationIndex()
    {
        return combinationIndex;   
    }


    void buildParameterCombinations()
    {
        if (cartesianProduct)
        {
            if (numberOfParameters == 0)
            {
                parameterCombinations = {};
            }
            else if (numberOfParameters == 1)
            {
                for (size_t i = 0; i < paramSpace[0].size(); ++i)
                {
                    std::vector<double> params = {paramSpace[0][i]};
                    parameterCombinations.push_back(params);
                }
            }
            else if (numberOfParameters == 2)
            {
                for (size_t i = 0; i < paramSpace[0].size(); ++i)
                {
                    for (size_t j = 0; j < paramSpace[1].size(); ++j)
                    {
                        std::vector<double> params = {paramSpace[0][i], paramSpace[1][j]};
                        parameterCombinations.push_back(params);
                    }
                }
            }
            else if (numberOfParameters == 3)
            {
                for (size_t i = 0; i < paramSpace[0].size(); ++i)
                {
                    for (size_t j = 0; j < paramSpace[1].size(); ++j)
                    {
                        for (size_t k = 0; k < paramSpace[2].size(); ++k)
                        {
                            std::vector<double> params = {paramSpace[0][i], paramSpace[1][j], paramSpace[2][k]};
                            parameterCombinations.push_back(params);
                        }
                    }
                }
            }
        }
        else
        {
            parameterCombinations = paramSpace;
            numberOfParameters = parameterCombinations[0].size();
            std::cout << "index " << index << std::endl;
            std::cout << "number of params " << numberOfParameters << std::endl; 
        }
    }


    void removeParams(const std::vector<std::vector<double>>& paramsToRemove)
    {
        parameterCombinations.erase(
        std::remove_if(parameterCombinations.begin(), parameterCombinations.end(),
            [&paramsToRemove](const std::vector<double>& param) {
                return std::find(paramsToRemove.begin(), paramsToRemove.end(), param) != paramsToRemove.end();
            }),
            parameterCombinations.end()
        );
    }


    double callFunction(double x) 
    {  
        std::vector<double> params;
        if (numberOfParameters != 0)
        {
            params = getParameterCombinationByIndex(combinationIndex);
        }
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
                return (*functionPtr)(x, params[0]);
            }
        }
        else if (index == 3 && numberOfParameters == 2)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, params[0], params[1]);
            }
        }
        else if (index == 4 && numberOfParameters == 3)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double& , const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, params[0], params[1], params[2]);
            }
        }
        else if (index == 0)
        {
            return callFunctionPython(x, params);
        }
        else
        {
            throw std::runtime_error("Number of parameters is not consistent with passed function");
        }
    }


    double callFunctionPython(double x, const std::vector<double>& params) 
    {
        auto function = std::get_if<py::function>(&functionVariant);
        if (numberOfParameters == 0) 
        {
            return (*function)(x).template cast<double>();  
        } 
        else if (numberOfParameters == 1)
        {
            double param1 = params[0];
            return (*function)(x, param1).template cast<double>();
        }
        else if (numberOfParameters == 2)
        {
            double param1 = params[0];
            double param2 = params[1];
            return (*function)(x, param1, param2).template cast<double>();
        }
        else if (numberOfParameters == 3)
        {
            double param1 = params[0];
            double param2 = params[1];
            double param3 = params[2];
            return (*function)(x, param1, param2, param3).template cast<double>();
        }
    }


    double getNumberOfParameters()
    {
        return numberOfParameters;
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    std::vector<std::vector<double>> getParameterCombinations() const
    {
        return parameterCombinations;
    }


    size_t getParameterCombinationsSize() const
    {
        return parameterCombinations.size();
    } 


    std::vector<double> getParameterCombinationByIndex(size_t index) const 
    {
        return parameterCombinations[index];
    }
};