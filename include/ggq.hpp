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
#include "optimizer.hpp"
#include "evaluator.hpp"
//#include "legendre_precompute.hpp"
using namespace std::numbers;


template<typename T>
class QuadratureRule
{
protected:
    size_t k;
    double lowerBound;
    double upperBound;
    double discretizerPrecision;
    double quadraturePrecision;
    T handler;
    matrix<double> invertedLegendreMatrix;
    std::vector<double> consolidatedEndpoints;
    //Legendre nodes and weights - rename
    std::vector<double> nodes;
    std::vector<double> weights;
    //values referes to values of the original function with each parameter - rename this
    std::vector<std::vector<double>> values;
    std::vector<std::vector<double>> compressedBasis;
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis; 
    std::vector<std::vector<double>> splitNodes;
    std::vector<double> chebyshevNodes; 
    std::vector<double> chebyshevWeights;
    std::vector<std::vector<std::vector<double>>> basisCoefficients;
    std::vector<double> basisIntegrals;
    std::vector<std::vector<double>> Jacobian; 


public:
    QuadratureRule(double lowerBoundInput, double upperBoundInput, T handler, double discretizerPrecisionInput = 1e-6, 
    double quadraturePrecisionInput = 1e-4, size_t kInput = 30) 
    : lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handler), discretizerPrecision(discretizerPrecisionInput),
    quadraturePrecision(quadraturePrecisionInput), k(kInput)
    {
        validatePrecisions(discretizerPrecision, quadraturePrecision);
    }


    ~QuadratureRule() 
    {

    }


    void calculateQuadratureNodes()
    {
        //Stage 1 of the paper 
        calculateConsolidatedEndpoints();
        sortConsolidatedEndpoints();
        removeDuplicateEndpoits();
        determineFinalNodes();
        determineValues();
    }


    void compressFunctionSpace() 
    {   
        //Stage 2 of the paper
        Compressor compressor(this, quadraturePrecision);
        compressedBasis = compressor.getCompressedBasis();
        chebyshevNodes = compressor.getChebyshevNodes();
        chebyshevWeights = compressor.getChebyshevWeights();
        splitNodesAtEndpoints(k);
    }


    void obtainBasisCoefficients()
    {
        splitCompressedBasisAtEndpoints();
        interpolateBasisFunctionsAtEachInterval();
    }


    void optimizeQuadrature()
    {
        //Stage 3 of the paper
        evaluateBasisIntegrals();
        Optimizer optimizer(chebyshevNodes, chebyshevWeights, basisCoefficients, splitCompressedBasis, basisIntegrals, consolidatedEndpoints, 
        splitNodes);
    }


    std::vector<std::vector<double>> evaluateAtChebyshevNodes()
    {
        size_t sizeOfParameterCombinations = handler.getParameterCombinationsSize();
        std::vector<std::vector<double>> valuesAll;

        for (size_t i = 0; i < sizeOfParameterCombinations; ++i)
        {   
            std::vector<double> values(chebyshevNodes.size());
            std::transform(chebyshevNodes.begin(), chebyshevNodes.end(), values.begin(), [this](double value){
                return this->handler.callFunction(value);
            });
            valuesAll.push_back(values);
            T::incrementCombinationIndex();
        }

        T::resetCombinationIndex();

        return valuesAll;
    }


    std::vector<double> evaluateIntegralsChebyshevNodes() 
    {
        std::vector<std::vector<double>> values = evaluateAtChebyshevNodes();

        std::vector<double> integrals(values.size());

        for (size_t i = 0; i < values.size(); ++i)
        {
            double integral;
            integral = innerProduct(values[i], chebyshevWeights);
            integrals[i] = integral;
        }

        return integrals;
    }


    std::vector<double> evaluateIntegrals()
    {
        std::vector<double> integrals(values.size());

        for (size_t i = 0; i < values.size(); ++i)
        {
            double integral;
            integral = innerProduct(values[i], weights);
            integrals[i] = integral;
        }

        return integrals;
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


    std::vector<std::vector<double>> getCompressedBasis() const
    {
        return compressedBasis;
    }


    std::vector<double> getChebyshevNodes() const
    {
        return chebyshevNodes;
    }


    std::vector<double> getChebyshevWeights() const
    {
        return chebyshevWeights;
    }


    std::vector<std::vector<double>> getSplitNodes() const
    {
        return splitNodes;
    }

    
    std::vector<std::vector<std::vector<double>>> getSplitCompressedBasis() const 
    {
        return splitCompressedBasis; 
    }


    std::vector<std::vector<std::vector<double>>> getBasisCoefficients() const
    {
        return basisCoefficients;
    }


    std::vector<double> getBasisFunctionInterval(const size_t& functionIndex = 0, const size_t& nodesIndex = 0)
    {
        std::vector<double> result;

        size_t indexStart = nodesIndex * splitNodes[nodesIndex].size();

        for (size_t i = indexStart; i < indexStart + splitNodes[nodesIndex].size(); ++i)
        {
            result.push_back(compressedBasis[functionIndex][i]);
        } 

        return result;
    }


    std::vector<double> getBasisIntegrals() const
    {
        return basisIntegrals; 
    } 


    std::vector<std::vector<double>> getJacobian() const
    {
        return Jacobian;
    }


protected:
    void evaluateBasisIntegrals()
    {
        basisIntegrals.resize(compressedBasis.size());

        for (size_t i = 0; i < compressedBasis.size(); ++i)
        {
            basisIntegrals[i] = innerProduct(compressedBasis[i], weights);
        }
    }


    void interpolateBasisFunctionsAtEachInterval()
    {
        basisCoefficients.resize(splitCompressedBasis.size());

        for (size_t i = 0; i < splitCompressedBasis.size(); ++i)
        {
            basisCoefficients[i].resize(splitCompressedBasis[i].size());

            for (size_t j = 0; j < splitCompressedBasis[i].size(); ++j)
            {
                Interpolator interpolator(k / 2, consolidatedEndpoints[j], consolidatedEndpoints[j+1], handler, splitCompressedBasis[i][j]);
                interpolator.interpolateFunction();

                basisCoefficients[i][j] = interpolator.getAlphaVector();
            }
        }
    }

    void splitCompressedBasisAtEndpoints() 
    {  
        splitCompressedBasis.resize(compressedBasis.size());

        for (size_t i = 0; i < compressedBasis.size(); ++i)
        {
            size_t numberOfIntervals = splitNodes.size(); 
            size_t valuesIndex = 0;

            for (size_t j = 0; j < numberOfIntervals; ++j)
            {
                size_t tempK = 0;
                std::vector<double> uValue; 

                while (tempK < k)
                {
                    uValue.push_back(compressedBasis[i][valuesIndex + tempK]);
                    tempK += 1;
                }

                valuesIndex += k;
                splitCompressedBasis[i].push_back(uValue);
            }
        }        
    }


    void splitNodesAtEndpoints(const size_t& k)
    {
        size_t valuesIndex = 0;

        for (size_t i = 0; i < consolidatedEndpoints.size() - 1; ++i)
        {
            std::vector<double> x = getNodesBetweenBounds(nodes, consolidatedEndpoints[i], consolidatedEndpoints[i+1]);

            splitNodes.push_back(x);
            
            valuesIndex += k;
        }
    }


    std::vector<double> getNodesBetweenBounds(const std::vector<double> x, const double& lowerBound, const double& upperBound)
    {
        std::vector<double> result;

        for (size_t i = 0; i < x.size(); ++i)
        {
            if (x[i] >= lowerBound && x[i] <= upperBound)
            {
                result.push_back(x[i]);
            }
        }

        return result;
    }


    //After a set of nodes are already stored, those need to be interpolated.
    //Where is the functionality of evaluating nodes and getting alphas? The only difference is that 
    //nodes are not dynamically determined but are taken as an input to Discretizer. 
    //What's the alogrithm? Input the endpoints check if they are already precise enough. Only if not start from the beginning.

    void calculateConsolidatedEndpoints()
    {
        size_t sizeOfParameterCombinations = handler.getParameterCombinationsSize();

        for (size_t i = 0; i < sizeOfParameterCombinations; ++i)
        {
            Discretizer<T> discretizer(k, discretizerPrecision, lowerBound, upperBound, handler);
            if (consolidatedEndpoints.size() != 0 && discretizer.evaluateStoppingConditionOfExistingEndpoints(consolidatedEndpoints))
            {
                std::cout << "Discretize function number " << i + 1 << std::endl;
                T::incrementCombinationIndex();
                continue;
            }
            else 
            {
                discretizer.discretize();
                std::vector<double> endpoints = discretizer.getFinalEndpoints();
                consolidatedEndpoints.insert(consolidatedEndpoints.end(), endpoints.begin(), endpoints.end());
                sortConsolidatedEndpoints();
                removeDuplicateEndpoits();
                std::cout << "Discretize function number " << i + 1 << std::endl;
                T::incrementCombinationIndex();
            }
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
        for (size_t i = 0; i < consolidatedEndpoints.size() - 1; ++i)
        {
            Interpolator divider(k / 2, consolidatedEndpoints[i], consolidatedEndpoints[i+1], handler);
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