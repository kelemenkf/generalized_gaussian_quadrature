#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include "utils.hpp"
#include "discretizer.hpp"

template <typename InputClass>
class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;

    using InputMethodType = double (InputClass::*)(const double&);
    InputClass* objectPtr;
    InputMethodType methodPtr;
    std::function<double(const double&)> method;

    using InputFunctionType = double(*)(const double&);
    InputFunctionType functionPtr;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType function = nullptr, 
    InputMethodType inputMethod = nullptr, InputClass* inputObject = nullptr) : lowerBound(lowerBoundInput), upperBound(upperBoundInput),
    functionPtr(function), methodPtr(inputMethod), objectPtr(inputObject)
    {
        validateFunctionExistence();
        if (objectPtr && methodPtr)
            method = std::bind(methodPtr, objectPtr, std::placeholders::_1);
    }

    ~QuadratureRule() 
    {

    }


    void discretizeFunctions()
    {
        int k = 30;
        double precision = 0.01;
        Discretizer discretizer(k, precision, lowerBound, upperBound, functionPtr, nullptr, nullptr);

        std::vector<double> endpoints;
        endpoints = discrteizer.discretizationRoutine();
    }


private: 
    void validateFunctionExistence()
    {
        if (functionPtr == nullptr && methodPtr == nullptr && objectPtr == nullptr)
        {
            throw std::invalid_argument("No function was supplied");
        }
    }


protected:
    double methodCaller(const double& value)
    {
        return method(value);
    }
};
#endif