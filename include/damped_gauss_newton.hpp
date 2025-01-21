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
    std::vector<double> x_c;
    MatrixXd initialJacobian; 
    unsigned int numberOfSteps; 
    double alpha; 
    double lambda; 
    double beta;


public:
    DGN()
    {};

    ~DGN() {};

private: 

};


#endif 