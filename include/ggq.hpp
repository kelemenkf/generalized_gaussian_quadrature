#ifndef GGQ_HPP
#define GGQ_HPP

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numbers>
#include "function_handler.hpp"
#include "discretizer.hpp"
#include "compressor.hpp"
using namespace std::numbers;


template<typename T>
class QuadratureRule
{
protected:
    double lowerBound;
    double upperBound;
    double discretizerPrecision;
    double quadraturePrecision;
    T handler;
    std::vector<double> consolidatedEndpoints;
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;
    std::vector<std::vector<double>> compressedBasis;


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, T handler, double discretizerPrecisionInput = 1e-6, 
    double quadraturePrecisionInput = 1e-4) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler), discretizerPrecision(discretizerPrecisionInput),
    quadraturePrecision(quadraturePrecisionInput)
    {
        validatePrecisions(discretizerPrecision, quadraturePrecision);
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
        determineValues();
    }


    void compressFunctionSpace() 
    {   
        Compressor(this, quadraturePrecision);
        std::cout << "It works" << std::endl;
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


    std::vector<std::vector<double>> getValues() const
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
        size_t sizeOfParameterCombinations = handler.getParameterCombinationsSize();

        for (size_t i = 0; i < sizeOfParameterCombinations; ++i)
        {
            Discretizer<T> discretizer(k, discretizerPrecision, lowerBound, upperBound, handler);
            std::vector<double> endpoints = discretizer.getFinalEndpoints();
            consolidatedEndpoints.insert(consolidatedEndpoints.end(), endpoints.begin(), endpoints.end());
            std::cout << "Discretize function number " << i + 1 << std::endl;
            T::incrementCombinationIndex();
        }

        T::resetCombinationIndex();
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
        size_t sizeOfParameterCombinations = handler.getParameterCombinationsSize();

        for (size_t j = 0; j < sizeOfParameterCombinations; ++j)
        {
            std::vector<double> functionValues(nodes.size());
            std::transform(nodes.begin(), nodes.end(), functionValues.begin(), [this](double value){
                return this->handler.callFunction(value);
            });
            values.push_back(functionValues);
            T::incrementCombinationIndex();
        }

        T::resetCombinationIndex();
    }


    void validatePrecisions(double discretizer, double quadrature)
    {
        if (quadrature / discretizer < 100)
        {
            throw std::invalid_argument("Quadrature precision has to be at least 100x of discretizer precision");
        }
    }
};

#endif