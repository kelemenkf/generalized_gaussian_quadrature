#include <iostream>
#include <functional>
#include <algorithm>
#include <cstdarg> 
#include <vector>

template<typename... Args>
class FunctionHandler
{
private:
    using InputFunction =  std::variant
    <
        std::function<double(const double&)>,
        std::function<double(const double&, const double&, const double&, const double&)>
    >;

    InputFunction functionVariant;
    std::vector<std::vector<double>> paramSpace; 
    size_t numberOfParameters;


public:
    FunctionHandler(InputFunction inputFunction, Args... args) : functionVariant(inputFunction) 
    {
        (paramSpace.push_back(args),...);
        numberOfParameters = sizeof...(args);
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }


    double callFunction(double x, double param1 = 0, double param2 = 0, double param3 = 0) 
    {
        if (numberOfParameters == 1) 
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&)>>(&functionVariant)) 
            {
                return (*functionPtr)(x);
            }
        } 
        else if (numberOfParameters == 3)
        {
            if (auto functionPtr = std::get_if<std::function<double(const double&, const double&, const double&, const double&)>>(&functionVariant))
            {
                return (*functionPtr)(x, param1, param2, param3);
            }
        }
    }


private:
};