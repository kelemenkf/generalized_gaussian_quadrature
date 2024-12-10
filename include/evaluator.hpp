#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP 
#include "utils.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

using namespace boost::math::tools;
using namespace boost::numeric::ublas;


class Evaluator
{
private: 
    const std::vector<double>& nodes;
    const std::vector<double>& values;
    const std::vector<double>& endpoints;
    std::vector<vector<double>> coefficients;


public: 
    Evaluator(const std::vector<double>& inputNodes, const std::vector<double>& inputValues, const std::vector<double>& inputEndpoints) :
    values(inputValues), nodes(inputNodes), endpoints(inputEndpoints)
    {
        //storeLagrangeCoefficientsOverAllIntervals();
    }

    ~Evaluator() {};


    std::vector<vector<double>> getCoefficients() const
    {
        return coefficients;
    }


    matrix<double> getLegendreMatirx() 
    {
        std::vector<double> x = divideNodesIntoIntervals(nodes, endpoints[0], endpoints[1]);
        matrix<double> L = calculateLegendrePolynomials(x);

        return L; 
    }

    
    double evaluate(double x, vector<double> coefficients)
    {
        double result = 0; 

        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            result += coefficients[i] * boost::math::legendre_p(i, x);
        }

        return result;
    }


private: 
    void storeLagrangeCoefficientsOverAllIntervals()
    {
        size_t valuesIndex = 0;

        for (size_t i = 0; i < endpoints.size() - 1; ++i)
        {
            std::vector<double> x = divideNodesIntoIntervals(nodes, endpoints[i], endpoints[i+1]);

            size_t intervalLength = x.size();

            std::vector<double> y;

            for (size_t j = valuesIndex; j < (valuesIndex + intervalLength); ++j)
            {
                y.push_back(values[j]);
            }

            matrix<double> L = calculateLegendrePolynomials(x);

            matrix<double> invertedL = invertMatrix(L);

            coefficients.push_back(calculateAlphaCoefficients(y, invertedL));

            valuesIndex += 30;
        }
    }


    std::vector<double> divideNodesIntoIntervals(const std::vector<double> x, const double& lowerBound, const double& upperBound)
    {
        std::vector<double> result;

        for (size_t i = 0; i < x.size(); ++i)
        {
            if (x[i] >= lowerBound && x[i] <= upperBound)
            {
                result.push_back(x[i]);
            }
        }

        return result;
    }



    polynomial<double> calculateLagrangeCoefficients(const std::vector<double>& x, const std::vector<double>& y)
    {
        size_t n = x.size();
        polynomial<double> result(0.0);

        for (size_t i = 0; i < n; ++i)
        {
            polynomial<double> term(1.0);
            double denom = 1.0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    term *= polynomial<double>({-x[j], 1.0});
                    denom *= (x[i] - x[j]);
                }
            }

            result += (term/denom) * y[i];
        }

        return result;
    }


    matrix<double> calculateLegendrePolynomials(const std::vector<double>& x)
    {
        matrix<double> legendreMatrix(x.size(), x.size());

        for (size_t i = 0; i < x.size(); ++i)
        {
            for (size_t j = 0; j < x.size(); ++j)
            {  
                if (j == 0)
                {
                    legendreMatrix(i, j) = 1;
                }
                else
                {
                    legendreMatrix(i, j) = boost::math::legendre_p(j, x[i]);
                }
            }
        }

        return legendreMatrix;
    }


    matrix<double> invertMatrix(matrix<double>& input) 
    {
        matrix<double> invertedLegendreMatrix(input.size1(), input.size2());

        matrix<double> copyOfLegendreMatrix(input);
        permutation_matrix<std::size_t> pm(copyOfLegendreMatrix.size1());

        int res = lu_factorize(copyOfLegendreMatrix, pm);
        if (res != 0) {
            throw std::invalid_argument("Matrix is singular");
        }

        invertedLegendreMatrix.assign(identity_matrix<double>(copyOfLegendreMatrix.size1()));
    
        lu_substitute(copyOfLegendreMatrix, pm, invertedLegendreMatrix);

        return invertedLegendreMatrix;
    }


    vector<double> calculateAlphaCoefficients(const std::vector<double>& P, const matrix<double>& invertedL)
    {
        vector<double> coefficients; 
        vector<double> lagrangeVectorUblas;
        lagrangeVectorUblas.resize(P.size());

        for (size_t i = 0; i < P.size(); ++i)
        {
            lagrangeVectorUblas(i) = P[i];
        }
        
        coefficients = prod(invertedL, lagrangeVectorUblas);

        return coefficients;
    }
};


#endif