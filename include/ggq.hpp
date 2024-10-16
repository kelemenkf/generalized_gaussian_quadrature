#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include "utils.hpp"
#include "function_handler.hpp"


class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;

    using InputFunctionType = std::function<double(const double&)>;
    InputFunctionType function;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType functionInput) : lowerBound(lowerBoundInput), upperBound(upperBoundInput),
    function(functionInput)
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
        if (!function) {
            throw std::runtime_error("Function is not initialized.");
        }
    }
};

#endif