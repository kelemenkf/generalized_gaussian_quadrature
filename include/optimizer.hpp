#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include "evaluator.hpp"
#include <Eigen/Dense>
#include <numeric>
#include "utils.hpp"
using namespace Eigen;


class Optimizer
{
private: 
    std::vector<double> chebyshevNodes;
    std::vector<double> chebyshevWeights;
    std::vector<std::vector<std::vector<double>>> basisCoefficients;
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis; 
    std::vector<double> basisIntegrals;
    MatrixXd Jacobian; 


public: 
    Optimizer(const std::vector<double>& inputChebyshevNodes, const std::vector<double>& inputChebyshevWeights, 
    const std::vector<std::vector<std::vector<double>>>& inputBasisCoefficients, const std::vector<std::vector<std::vector<double>>>& 
    inputSplitCompressedBasis, const std::vector<double>& inputBasisIntegrals)
    : chebyshevNodes(inputChebyshevNodes), chebyshevWeights(inputChebyshevWeights), basisCoefficients(inputBasisCoefficients), 
    splitCompressedBasis(inputSplitCompressedBasis), basisIntegrals(inputBasisIntegrals)
    { 
        validateSystem();
    };

    ~Optimizer(){};


protected: 
    void formJacobian()
    {
        Jacobian.resize(splitCompressedBasis.size(), 2 * chebyshevNodes.size());

        std::cout << Jacobian.rows() << " " << Jacobian.cols() << std::endl;
    }

private: 
    void validateSystem()
    {
        int m = basisIntegrals.size(); 
        int n = chebyshevNodes.size(); 

        if (m <= n)
        {
            throw std::invalid_argument("The system is not overdetermined.");
        }
    }
};


#endif