#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include "evaluator.hpp"
#include "damped_gauss_newton.hpp"
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
    std::vector<std::pair<double, int>> sortedStepDirectionNorms; 


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
        calculateStepDirectionNorms();
        reorderNodesBasedOnNorms();
        formInitialJacobianDGN(sortedStepDirectionNorms[0].second, chebyshevNodes.size()); 
        formX();
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


    std::vector<double> getStepDirectionNorms() const 
    {
        return stepDirectionNorms; 
    }


    std::vector<std::pair<double, int>> getSortedStepDirectionNorms() const 
    {
        return sortedStepDirectionNorms;
    }


protected:
    void dampedGaussNewton(int removedNode, int numberOfSteps)
    {
        std::cout << "Performing " << numberOfSteps << " number of damped Gauss-Newton steps." << std::endl;

        for (size_t i = 0; i < sortedStepDirectionNorms.size(); ++i)
        {
            //Calls DGN with all combinations of nodes removed in order of inreasing eta        
            //First take current nodes, remove the one on the current iteration then apply the calculated delta x_k. This is the frist step
            //and input for DGN. 
            VectorXd x_initial = stepDirections[sortedStepDirectionNorms[i].second];
            VectorXd x_1 = firstStep(x_initial, i, x_initial.size());
            //DGN();
        }
    }


    VectorXd firstStep(const VectorXd& x_initial, int i, int n)
    {
        VectorXd x = formX();   
        return removeNode(x, i, n) - removeNode(x_initial, i, n);
    }

    
    double calculatePrecision()
    {
        //takes a quadrature as an input and cacluates precision
        //for the original input this should already be known?  
    }


    void formInitialJacobianDGN(int columnToRemove, int n) 
    {
        MatrixXd JacobianDGN = removeColumnFromEigenMatrix(Jacobian, columnToRemove, n);

        std::cout << "System was reduced to a " << JacobianDGN.rows() << " by " << JacobianDGN.cols() << " system." << std::endl;
    }


    void reorderNodesBasedOnNorms()
    {
        for (size_t i = 0; i < stepDirectionNorms.size(); ++i)
        {
            sortedStepDirectionNorms.push_back({stepDirectionNorms[i], i});
        }
        
        std::sort(sortedStepDirectionNorms.begin(), sortedStepDirectionNorms.end());
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
        size_t n = chebyshevNodes.size(); 
        stepDirections.resize(n);

        for (size_t k = 0; k < chebyshevNodes.size(); ++k)
        {
            MatrixXd A_k = shermanMorrisonWoodburry(A, k, n);

            MatrixXd J_k = updateJacobian(Jacobian, k, n);

            VectorXd stepDirection = J_k.transpose() * A_k * eigenBasisIntegrals;

            stepDirections[k] = stepDirection;
        }
    }


    MatrixXd shermanMorrisonWoodburry(const MatrixXd& input, int k, int j = 0)
    {
		//this may be not good 
		
		MatrixXd A_prime = input;
	
		VectorXd u_k = Jacobian.col(k); 

		double scalar1 = 1.0 + u_k.transpose() * input * u_k;

		VectorXd Au_k = input * u_k; 

		A_prime -= (Au_k * Au_k.transpose()) / scalar1;

    	VectorXd u_kj = Jacobian.col(k + j);

		double scalar2 = 1.0 + u_kj.transpose() * A_prime * u_kj;

		VectorXd A_prime_u_kj = A_prime * u_kj;

		A_prime -= (A_prime_u_kj * A_prime_u_kj.transpose()) / scalar2;

		return A_prime;
    }


    MatrixXd updateJacobian(const MatrixXd& Jacobian, int k, int n)
    {
        VectorXd v(Jacobian.cols()); 
        v.setZero();
        v[k] = 1;

        VectorXd column = Jacobian.col(k);

        VectorXd v2(Jacobian.cols());
        v2.setZero();
        v2[k+n] = 1;

        VectorXd column2 = Jacobian.col(k + n);

        return Jacobian - (column * v.transpose()) - (column2 * v2.transpose());
    }


    void formA()
    {
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
                    Evaluator evaluator(alphas, lowerBound, upperBound, splitNodes[interval]);
                    double derivativeAtChebyshevNode = evaluator.evaluateFirstDerivative(((element->second)[i]).first);
                    double u = evaluator.evaluateUnreversed((element->second)[i].first);
                    Jacobian(function, nodeIndex) = derivativeAtChebyshevNode * (element->second)[i].second;  
                    Jacobian(function, nodeIndex + numberOfNodes) = u; 
                    nodeIndex += 1;
                }
            }
        }

        int m = Jacobian.rows(); 
        int n = Jacobian.cols(); 

        if (m < n)
        {
            std::cout << "System is underdetermined of size " << m << " by " << n << std::endl;
        }
        else if (m == n)
        {
            std::cout << "System is " << m << " by " << n << std::endl;
        }
    }


private:   
    VectorXd removeNode(VectorXd input, int node, int n)
    {
        std::cout << "Size of input vector before removal " << input.size() << std::endl;
        std::vector<double> temp = convertEigenVectorToStd<double>(input);

        std::vector<int> indicesToRemove = {node + n};
        indicesToRemove.push_back(node);

        for (size_t index: indicesToRemove)
        {
            if (index < temp.size())
            {
                temp.erase(temp.begin() + index);
            }
        }

        VectorXd result = convertStdVectorToEigen(temp); 

        std::cout << "Size of vector after removal " << result.size() << std::endl;

        return result;
    }

 
    VectorXd formX()
    {
        int n = chebyshevNodes.size();
        VectorXd x(2 * n);

        for (size_t node = 0; node < chebyshevNodes.size(); ++node)
        {
            x[node] = chebyshevNodes[node];
            x[node + n] = chebyshevWeights[node];
        }

        std::cout << "x for DGN formed at size " << x.size() << std::endl;
        std::cout << x << std::endl;

        return x;
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