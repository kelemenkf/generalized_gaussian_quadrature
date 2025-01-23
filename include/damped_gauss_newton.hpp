#ifndef DAMPED_GAUSS_NEWTON_HPP
#define DAMPED_GAUSS_NEWTON_HPP

#include "utils.hpp"
#include "evaluator.hpp"
#include <Eigen/Dense>
using namespace Eigen; 

class DGN
{
private: 
    std::vector<std::vector<std::vector<double>>> basisCoefficients;
    std::vector<double> endpoints;
    VectorXd R; 
    VectorXd x;
    MatrixXd Jacobian; 
    unsigned int numberOfSteps; 
    //is alpha an input? 
    double alpha; 
    double lambda; 
    double beta;
        


public:
    DGN(std::vector<std::vector<std::vector<double>>> inputBasisCoefficients, std::vector<double> inputEndpoints, VectorXd
    inputInitialX, MatrixXd inputInitialJacobian, unsigned int inputNumberOfSteps, double inputAlpha, double inputLambda = 10e-4, double inputBeta = 0.9)
    : basisCoefficients(inputBasisCoefficients), endpoints(inputEndpoints), x(inputInitialX), Jacobian(inputInitialJacobian), 
    numberOfSteps(inputNumberOfSteps), alpha(inputAlpha), lambda(inputLambda), beta(inputBeta)
    {};

    ~DGN() {};


private:
    void iterateOverNumberOfSteps()
    {
        unsigned int k = 0; 
        while (k < numberOfSteps)
        {
            updateJacobianAndR();
            //iterate();
            std::cout << "Damped Gauss Newton iteration number " << k + 1 << std::endl;
            ++k;
        }
    }


    void updateJacobianAndR()
    {
        for (size_t u = 0; u < basisCoefficients.size(); ++u)
        {
            double R_temp = 0; 
            for (size_t endpoint = 0; endpoint < endpoints.size(); ++endpoint)
            {
                int n = x.size() / 2; 
                for (int i = 0; i < n; ++i)
                {
                    if (x[i] > endpoints[endpoint] && x[i] < endpoints[endpoint+1])
                    {
                        Evaluator evaluator(basisCoefficients[u][endpoint], endpoints[endpoint], endpoints[endpoint+1], {});
                        double u_x = evaluator.evaluateUnreversed(x[i]); 
                        double u_x_prime = evaluator.evaluateFirstDerivative(x[i]);
                        double w = x[i + n];
                        R_temp += u*w;
                        Jacobian(u, i) = u_x_prime * w; 
                        Jacobian(u, i + n) = u_x;
                    }
                }
            }
            R[u] = R_temp;
        }
    }


    void iterate()
    {
        x = x - Jacobian.transpose() * (Jacobian * Jacobian.transpose() + lambda * MatrixXd::Identity(Jacobian.rows(), Jacobian.rows())) * R;
    }


    void checkWolfeConditons()
    {
    }
};


#endif 