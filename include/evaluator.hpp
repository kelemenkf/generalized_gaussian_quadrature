#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP 
#include "utils.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::math::tools;
using namespace boost::numeric::ublas;


class Evaluator
{
private: 
    //turn nodes into a map?
    std::vector<double> inputNodes; 
    std::vector<double> reversedNodes;
    std::vector<double> coefficients;
    double lowerBound; 
    double upperBound;
    std::vector<double> output; 


public:
    Evaluator(std::vector<double> inputValue, const std::vector<double>& coefficientsInput, double inputLowerBound, 
    double inputUpperBound) : inputNodes(inputValue), coefficients(coefficientsInput), lowerBound(inputLowerBound), upperBound(inputUpperBound) 
    {
        reverseNodes();
        displayVector(reversedNodes);
        evaluateAll();
    };


    ~Evaluator() {};


    std::vector<double> getOutput() const
    {
        return output; 
    }


    std::vector<double> getReversedNodes() const
    {
        return reversedNodes;
    }


    std::vector<double> getInputNodes() const
    {
        return inputNodes;
    }


private: 
    void evaluateAll()
    {
        output.resize(reversedNodes.size());

        std::transform(reversedNodes.begin(), reversedNodes.end(), output.begin(), [this](double value){
            return this->evaluate(value);
        });
    }


    double evaluate(double x)
    {
        double result = 0;

       for (size_t i = 0; i < coefficients.size(); ++i)
       {
            result += coefficients[i] * transformNode(boost::math::legendre_p(i, x));
       } 

       return result;
    }


    void reverseNodes()
    {
        reversedNodes.resize(inputNodes.size());

        std::transform(inputNodes.begin(), inputNodes.end(), reversedNodes.begin(), [this](double value){
            return this->reverseNode(value);
        });
    }


    inline double transformNode(double value) const
    {
        return ((this->upperBound - this->lowerBound) * value + (this->upperBound + this->lowerBound)) / 2;
    }


    inline double reverseNode(double value) const 
    {
        return (2 * value - (this->upperBound + this->lowerBound)) / (this->upperBound - this->lowerBound);
    }
};

#endif
