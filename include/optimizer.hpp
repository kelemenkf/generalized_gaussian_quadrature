#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include "evaluator.hpp"
#include <Eigen/Dense>
#include <numeric>
#include "utils.hpp"
#include <map>
using namespace Eigen;


class Optimizer
{
private: 
    std::vector<double> chebyshevNodes;
    std::vector<double> chebyshevWeights;
    std::vector<std::vector<std::vector<double>>> basisCoefficients;
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis; 
    std::vector<double> basisIntegrals;
    std::vector<double> endpoints;
    std::vector<std::vector<double>> splitNodes;
    std::map<int, std::vector<std::pair<double, double>>> intervalChebyshevNodesMap;
    MatrixXd Jacobian; 


public: 
    Optimizer(const std::vector<double>& inputChebyshevNodes, const std::vector<double>& inputChebyshevWeights, 
    const std::vector<std::vector<std::vector<double>>>& inputBasisCoefficients, const std::vector<std::vector<std::vector<double>>>& 
    inputSplitCompressedBasis, const std::vector<double>& inputBasisIntegrals, const std::vector<double>& inputEndpoints, 
    const std::vector<std::vector<double>>& inputSplitNodes)
    : chebyshevNodes(inputChebyshevNodes), chebyshevWeights(inputChebyshevWeights), basisCoefficients(inputBasisCoefficients), 
    splitCompressedBasis(inputSplitCompressedBasis), basisIntegrals(inputBasisIntegrals), endpoints(inputEndpoints), splitNodes(inputSplitNodes)
    { 
        assignChebyshevNodesToInterval();
        formJacobian();
    };

    ~Optimizer(){};


    std::map<int, std::vector<std::pair<double, double>>> getIntervalChebyshevNodesMap() const
    {
        return intervalChebyshevNodesMap;
    }


    std::vector<std::vector<double>> getJacobian()
    {
        std::vector<std::vector<double>> jacobian(Jacobian.rows());

        for (size_t row = 0; row < Jacobian.rows(); ++row)
        {
            jacobian[row].resize(Jacobian.cols());

            for (size_t col = 0; col < Jacobian.cols(); ++col)
            {
                jacobian[row][col] = Jacobian(row, col);
            }
        }

        return jacobian;
    }


protected: 
    void formJacobian()
    {
        Jacobian.resize(splitCompressedBasis.size(), 2 * chebyshevNodes.size());

        for (size_t function = 0; function < splitCompressedBasis.size(); ++function)
        {
            int nodeIndex = 0;
            int numberOfNodes = chebyshevNodes.size();

            for (auto element = intervalChebyshevNodesMap.begin(); element != intervalChebyshevNodesMap.end(); 
            ++element)
            {
                for(size_t i = 0; i < (element->second).size(); ++i)
                {
                    int interval = element->first;
                    int lowerBound = endpoints[interval]; 
                    int upperBound = endpoints[interval + 1]; 
                    std::cout << "Interval " << interval << " lower bound " << lowerBound << " upper bound " << upperBound << std::endl;
                    std::vector<double> alphas = basisCoefficients[function][interval];
                    displayVector(alphas);
                    Evaluator evaluator(splitNodes[interval], alphas, lowerBound, upperBound);
                    std::cout << (element->second)[i].first << " " << (element->second)[i].second << std::endl;
                    double derivativeAtChebyshevNode = evaluator.evaluateFirstDerivative(((element->second)[i]).first);
                    double u = evaluator.evaluateUnreversed((element->second)[i].first);
                    Jacobian(function, nodeIndex) = derivativeAtChebyshevNode * (element->second)[i].second;  
                    Jacobian(function, nodeIndex + numberOfNodes) = u; 
                    nodeIndex += 1;
                }
            }

        }

        std::cout << Jacobian << std::endl;
    }


private: 
    void sortChebyshevNodes()
    {
        std::sort(chebyshevNodes.begin(), chebyshevNodes.end());
    }


    void assignChebyshevNodesToInterval()
    {
       for (int interval = 0; interval < splitNodes.size(); ++interval)
       {
            std::vector<std::pair<double, double>> nodes;

            for (size_t node = 0; node < chebyshevNodes.size(); ++node)
            {
                if (chebyshevNodes[node] > endpoints[interval] && chebyshevNodes[node] < endpoints[interval + 1])
                {
                    nodes.push_back({chebyshevNodes[node], chebyshevWeights[node]});
                } 
            }

            intervalChebyshevNodesMap[interval] = nodes;
       } 
    }
};


#endif