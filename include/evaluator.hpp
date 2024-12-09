#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP 
#include "utils.hpp"


class Evaluator
{
private:
    std::vector<std::vector<double>> interpolationCoefficients;
    std::vector<std::vector<double>> nodes;


public: 
    Evaluator(std::vector<std::vector<double>> inputNodes, std::vector<std::vector<double>> inputInterpolationCoefficients) :
    interpolationCoefficients(inputInterpolationCoefficients), nodes(inputNodes)
    {}

    ~Evaluator() {};

private: 


};


#endif
