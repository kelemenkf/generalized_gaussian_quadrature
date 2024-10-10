#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>


template <typename InputClass>
class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;

    using InputMethodType = double (InputClass::*)(const double&);
    InputClass* objectPtr;
    InputMethodType methodPtr;

    using InputFunctionType = double(*)(const double&);
    InputFunctionType functionPtr;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType function = nullptr, 
    InputMethodType method = nullptr, InputClass* inputObject = nullptr) : lowerBound(lowerBoundInput), upperBound(upperBoundInput),
    functionPtr(function), methodPtr(method), objectPtr(inputObject)
    {
        validateFunctionExistence();
    };

    ~QuadratureRule() 
    {

    };


private: 
    void validateFunctionExistence()
    {
        if (functionPtr == nullptr && methodPtr == nullptr && objectPtr == nullptr)
        {
            throw std::invalid_argument("No function was supplied");
        }
    };
};

#endif