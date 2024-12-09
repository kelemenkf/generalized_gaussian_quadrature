#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP 
#include "utils.hpp"
#include <boost/math/tools/polynomial.hpp>

using namespace boost::math::tools;


class Evaluator
{
private:    
    const std::vector<double>& nodes;
    const std::vector<double>& values;
    polynomial<double> coefficients;


public: 
    Evaluator(const std::vector<double>& inputNodes, const std::vector<double>& inputValues) :
    values(inputValues), nodes(inputNodes)
    {
        calculateLagrangeCoefficients();
    }

    ~Evaluator() {};


    polynomial<double> getCoefficients() const
    {
        return coefficients;
    }


private: 
    void calculateLagrangeCoefficients()
    {
        size_t n = nodes.size();
        coefficients = polynomial<double>(0.0);

        displayVector(values);
        displayVector(nodes);

        for (size_t i = 0; i < n; ++i)
        {
            polynomial<double> term(1.0);
            double denom = 1.0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    term *= polynomial<double>({-nodes[j], 1.0});
                    denom *= (nodes[i] - nodes[j]);
                }
            }

            coefficients += (term/denom) * values[i];
        }
    }
};


#endif
