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
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, T handler) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler)
    {

    }

    ~QuadratureRule() 
    {

    }


    void calculateQuadratureNodes()
    {
        calculateConsolidatedEndpoints();
        sortConsolidatedEndpoints();
        removeDuplicateEndpoits();
        determineFinalNodes();
    }


    double getLowerBound() const 
    {
        return lowerBound;
    }


    std::vector<double> getConsolidatedEndpoints() const
    {
        return consolidatedEndpoints;
    }


    std::vector<double> getNodes() const
    {
        return nodes;
    }


    std::vector<double> getValues() const
    {
        return values;
    }


    std::vector<double> getWeights() const
    {
        return weights;
    }


protected:
    void calculateConsolidatedEndpoints()
    {
        int k = 30;
        double precision = 1e-5;
        size_t sizeOfParameterCombinations = handler.getParameterCombinationsSize();

        for (size_t i = 0; i < sizeOfParameterCombinations; ++i)
        {
            Discretizer<T> discretizer(k, precision, lowerBound, upperBound, handler);
            std::vector<double> endpoints = discretizer.getFinalEndpoints();
            consolidatedEndpoints.insert(consolidatedEndpoints.end(), endpoints.begin(), endpoints.end());
            std::cout << "Discretize function number " << i + 1 << std::endl;
        }
    }


    void sortConsolidatedEndpoints()
    {
        std::sort(consolidatedEndpoints.begin(), consolidatedEndpoints.end());
    }


    void removeDuplicateEndpoits()
    {
        auto last = std::unique(consolidatedEndpoints.begin(), consolidatedEndpoints.end());

        consolidatedEndpoints.erase(last, consolidatedEndpoints.end());
    }


    void determineFinalNodes()
    {
        int k = 30;
        for (size_t i = 0; i < consolidatedEndpoints.size() - 1; ++i)
        {
            IntervalDivider divider(k / 2, consolidatedEndpoints[i], consolidatedEndpoints[i+1], handler);
            divider.calculateLegendreNodes();
            std::vector<double> transformedNodes = divider.getTransformedMesh();
            nodes.insert(nodes.end(), transformedNodes.begin(), transformedNodes.end());
            std::vector<double> intervalWeights = divider.getQuadratureWeights();
            weights.insert(weights.end(), intervalWeights.begin(), intervalWeights.end());
        }
    }


    void determineValues()
    {
        
    }
};

#endif