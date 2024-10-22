#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numbers>
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
    std::vector<double> consolidatedEndpoints;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, T handler) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler)
    {

    }

    ~QuadratureRule() 
    {

    }

    void calculateConsolidatedEndpoints()
    {
        int k = 30;
        double precision = 1e-5;
        size_t sizeOfParameterCombinations = T::getParameterCombinationsSize();

        for (size_t i = 0; i < sizeOfParameterCombinations; ++i)
        {
            std::cout << T::getParameterCombinationsIndex() << std::endl;
            Discretizer<T> discretizer(k, precision, lowerBound, upperBound, handler);
            std::vector<double> endpoints = discretizer.getFinalEndpoints();
            consolidatedEndpoints.resize(consolidatedEndpoints.size() + endpoints.size());
            consolidatedEndpoints.insert(consolidatedEndpoints.end(), endpoints.begin(), endpoints.end());
            T::incrementCombinationIndex();
        }
    }


    double getLowerBound() const 
    {
        return lowerBound;
    }


    // std::vector<double> getConsolidatedEndpoints() const
    // {
    //     return consolidatedEndpoints;
    // }
};

#endif