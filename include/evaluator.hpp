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
    Evaluator(const std::vector<double>& coefficientsInput, double inputLowerBound, double inputUpperBound, std::vector<double> inputValues) : 
    inputNodes(inputValues), coefficients(coefficientsInput), lowerBound(inputLowerBound), upperBound(inputUpperBound) 
    {
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


    double evaluateFirstDerivative(const double& x)
    {
        double result = 0; 

        std::cout << x << " " << reverseNode(x) << std::endl;

        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            result += coefficients[i] * boost::math::legendre_p_prime(i, reverseNode(x));
        }

        return result;
    }


    double evaluateUnreversed(const double& x)
    {
        double result = 0; 

        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            double res = boost::math::legendre_p(i, reverseNode(x));
            result += coefficients[i] * res;
        }

        return result;
    }


    void evaluateInput()
    {
        reverseNodes();
        evaluateAll();
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
            double res = boost::math::legendre_p(i, x); 
            result += coefficients[i] * res;
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
