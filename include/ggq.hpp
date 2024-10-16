#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include "utils.hpp"
#include "function_handler.hpp"


template<typename T>
class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;
    T handler;

public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, const T& handler) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler)
    {
        
    }

    ~QuadratureRule() 
    {

    }

    void discretizeFunctions()
    {
    }
};

#endif