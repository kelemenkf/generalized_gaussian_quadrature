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
    VectorXd r; 
    VectorXd x;
    MatrixXd Jacobian; 
    unsigned int numberOfSteps; 
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
            evaluateRAtNewNodes();
            modifyJacobian();
            iterate();
        }
    }


    void modifyJacobian()
    {

    }


    void evaluateRAtNewNodes()
    {
        Evaluator(basisCoefficients);
    }


    void iterate()
    {
        x = x - Jacobian.transpose() * (Jacobian * Jacobian.transpose() + MatrixXd::Identity(Jacobian.rows(), Jacobian.rows())) * r;
    }


    void checkWolfeConditons()
    {
        //Does it need evaluation, if so 
    }
};


#endif 