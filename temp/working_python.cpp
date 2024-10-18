#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdarg> 
#include <vector>

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


public:
    FunctionHandler(InputFunction inputFunction, Parameter... parameters) : functionVariant(inputFunction) 
    {
        (paramSpace.push_back(parameters),...);
        numberOfParameters = sizeof...(parameters);
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    double callFunction(double x, double param1 = 0, double param2 = 0, double param3 = 0) 
    {
        size_t index = functionVariant.index(); 
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
        else if (index == 3 & numberOfParameters == 3)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2, param3);
            }
        }
        else
        {
            throw std::runtime_error("Number of parameters is not consistent with passed function");
        }
    }
};
