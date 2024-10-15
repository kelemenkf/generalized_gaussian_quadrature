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

public:
    FunctionHandler(Args... args) 
    {
        (paramSpace.push_back(args),...);
    }


    std::vector<std::vector<double>> getParamSpace() const
    {
        return paramSpace;
    }

private:
};