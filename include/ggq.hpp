#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>

template <typename InputClass>
class QuadratureRule
{
private:
    double lowerBound;
    double upperBound;

    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    InputClass* objectPtr;
    InputMethodType methodPtr;

    using InputFunctionType = double(*)(const std::vector<double>&);
    InputFunctionType functionPtr;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType function = nullptr, 
    InputMethodType mehtod = nullptr, InputClass* inputObject = nullptr)
    {
        validateFunctionExistence();
    };

    ~QuadratureRule() 
    {

    };


private: 
    static double validateLowerBound(double input);

    static double validateUpperBound(double input);

    void validateFunctionExistence()
    {
        if (!((functionPtr) || (methodPtr && objectPtr))) 
        {
            throw std::invalid_argument("No function was supplied");
        }
    };
};

#endif