#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numbers>
#include "utils.hpp"
#include "function_handler.hpp"
#include "discretizer.hpp"
using namespace std::numbers;


template<typename T>
class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;
    T handler;

public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, T handler) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler)
    {

    }

    ~QuadratureRule() 
    {

    }

    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> discretizeFunctions()
    {
        int k = 30;
        double precision = 1e-5;

        Discretizer<T> discretizer(k, precision, lowerBound, upperBound, handler);

        return discretizer.determineFinalNodes();
    }


    double getLowerBound() const 
    {
        return lowerBound;
    }
};

#endif