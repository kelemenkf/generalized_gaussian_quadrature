#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include "utils.hpp"

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
    QuadratureRule(double lowerBoundInput = -1, double upperBoundInput = 1, InputFunctionType function = nullptr, 
    InputMethodType inputMethod = nullptr, InputClass* inputObject = nullptr) : lowerBound(lowerBoundInput), upperBound(upperBoundInput),
    functionPtr(function), methodPtr(inputMethod), objectPtr(inputObject)
    {
        validateFunctionExistence();
    }

    ~QuadratureRule() 
    {

    }


    void discretizeFunctions()
    {
    }


private: 
    void validateFunctionExistence()
    {
        if (functionPtr == nullptr && methodPtr == nullptr && objectPtr == nullptr)
        {
            throw std::invalid_argument("No function was supplied");
        }
    }
};

#endif