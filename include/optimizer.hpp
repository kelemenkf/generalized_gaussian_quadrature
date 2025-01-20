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
    VectorXd eigenBasisIntegrals;
    std::vector<double> endpoints;
    std::vector<std::vector<double>> splitNodes;
    std::map<int, std::vector<std::pair<double, double>>> intervalChebyshevNodesMap;
    MatrixXd Jacobian; 
    MatrixXd A;
    std::vector<VectorXd> stepDirections; 
    std::vector<double> stepDirectionNorms;


public: 
    Optimizer(const std::vector<double>& inputChebyshevNodes, const std::vector<double>& inputChebyshevWeights, 
    const std::vector<std::vector<std::vector<double>>>& inputBasisCoefficients, const std::vector<std::vector<std::vector<double>>>& 
    inputSplitCompressedBasis, const std::vector<double>& inputBasisIntegrals, const std::vector<double>& inputEndpoints, 
    const std::vector<std::vector<double>>& inputSplitNodes)
    : chebyshevNodes(inputChebyshevNodes), chebyshevWeights(inputChebyshevWeights), basisCoefficients(inputBasisCoefficients), 
    splitCompressedBasis(inputSplitCompressedBasis), basisIntegrals(inputBasisIntegrals), endpoints(inputEndpoints), splitNodes(inputSplitNodes)
    { 
        assignChebyshevNodesToInterval();
        transformIntegrals();
        formJacobian();
        formA();
        calculateStepDirections();
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


    std::vector<VectorXd> getStepDirections() const 
    {
        return stepDirections;
    }


protected:
    void reorderNodesBasedOnNorms()
    {
        
    }


    void calculateStepDirectionNorms()
    {
        stepDirectionNorms.resize(stepDirections.size());

        for (size_t i = 0; i < stepDirections.size(); ++i)
        {
            stepDirectionNorms[i] = stepDirections[i].norm();
        }
    }


    void transformIntegrals()
    {
        eigenBasisIntegrals.resize(basisIntegrals.size());

        for (size_t i = 0; i < basisIntegrals.size(); ++i)
        {
            eigenBasisIntegrals(i) = basisIntegrals[i];
        }

    }


    void calculateStepDirections()
    {
        size_t j = chebyshevNodes.size(); 
        stepDirections.resize(j);

        std::cout << "a" << std::endl;

        for (size_t k = 0; k < chebyshevNodes.size(); ++k)
        {
            MatrixXd A_k = shermanMorrisonWoodburry(A, k, j);

            std::cout << "A_k size " << A_k.rows() << " " << A_k.cols() << std::endl;

            std::cout << "Jacobian size " << Jacobian.rows() << " " << Jacobian.cols() << std::endl;

            std::cout << "r size" << eigenBasisIntegrals.rows() << " " << eigenBasisIntegrals.cols() << std::endl;

            VectorXd stepDirection = A_k * (Jacobian.transpose() * eigenBasisIntegrals);

            stepDirections[k] = stepDirection;
        }
    }


    MatrixXd shermanMorrisonWoodburry(const MatrixXd& input, int k, int j = 0)
    {
		//this may be not good 
		
		MatrixXd A_prime = input;
	
		VectorXd u_k = Jacobian.col(k); 

        std::cout << "h" << std::endl;

		double scalar1 = 1.0 + u_k.transpose() * input * u_k;

        std::cout << "g" << std::endl;

		VectorXd Au_k = input * u_k; 

        std::cout << "f" << std::endl;

		A_prime -= (Au_k * Au_k.transpose()) / scalar1;

    	VectorXd u_kj = Jacobian.col(k + j);

        std::cout << "e" << std::endl;

		double scalar2 = 1.0 + u_kj.transpose() * A_prime * u_kj;

        std::cout << "d" << std::endl; 

		VectorXd A_prime_u_kj = A_prime * u_kj;

        std::cout << "b" << std::endl;

		A_prime -= (A_prime_u_kj * A_prime_u_kj.transpose()) / scalar2;

        std::cout << "c" << std::endl;

		return A_prime;
    }


    void formA()
    {
        std::cout << "ez meg jo" << std::endl;
        A = (Jacobian * Jacobian.transpose()).inverse();
    }


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
                    std::vector<double> alphas = basisCoefficients[function][interval];
                    Evaluator evaluator(splitNodes[interval], alphas, lowerBound, upperBound);
                    double derivativeAtChebyshevNode = evaluator.evaluateFirstDerivative(((element->second)[i]).first);
                    double u = evaluator.evaluateUnreversed((element->second)[i].first);
                    Jacobian(function, nodeIndex) = derivativeAtChebyshevNode * (element->second)[i].second;  
                    Jacobian(function, nodeIndex + numberOfNodes) = u; 
                    nodeIndex += 1;
                }
            }

        }
    }


private: 


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