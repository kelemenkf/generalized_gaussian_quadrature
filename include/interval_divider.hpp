#ifndef INTERVAL_DIVIDER_HPP
#define INTERVAL_DIVIDER_HPP

#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <cmath>
#include "utils.hpp"
using namespace boost::numeric::ublas;


template<typename T>
class IntervalDivider
{
private: 
    int k;
    std::vector<double> legendreMesh;
    std::vector<double> transformedMesh;
    std::vector<double> quadratureWeights;
    matrix<double> legendreMatrix;
    matrix<double> invertedLegendreMatrix;
    vector<double> alphaVector;
    std::vector<double> lagrangeVector;
    double measure;
    double lowerBound;
    double upperBound;
    T handler;


public: 
    IntervalDivider(int kInput, double lowerBoundInput, double upperBoundInput, const T& handlerInput, const std::vector<double>& 
    lagrangeInput = {})  
    : k(kInput), lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handlerInput), lagrangeVector(lagrangeInput)
    {};


    ~IntervalDivider() {};


    void processInterval()
    {
        calculateMesh();
        transformMesh();
        calculateLegendrePolynomials(legendreMatrix);
        invertMatrix();
        evaluateFunctionOnTransformedMesh();
        calculateAlphaCoefficients();
        calculateSquaredAlphas();
    }


    void calculateLegendreNodes()
    {
        calculateMesh();
        transformMesh();
        calculateWeights();
    }


    void interpolateFunction()
    {
        calculateMesh();
        transformMesh();
        calculateLegendrePolynomials(legendreMatrix);
        invertMatrix();
        calculateAlphaCoefficients();
    }


    double evaluate(double x, vector<double> coefficients)
    {
        double result = 0; 

        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            std::cout << coefficients[i] << " " << result << " " << boost::math::legendre_p(i, x) << std::endl;
            result += coefficients[i] * boost::math::legendre_p(i, x);
        }

        return result;
    }


    double getMeasure() const
    {
        return measure;
    }

    
    std::vector<double> getLegendreMesh() const
    {
        return legendreMesh;
    }

    
    std::vector<double> getTransformedMesh() const
    {
        return transformedMesh;
    }


    std::vector<double> getQuadratureWeights() const 
    {
        return quadratureWeights;
    }


    matrix<double> getLegendreMatrix() const
    {
        return legendreMatrix;
    }


    matrix<double> getInvertedLegendreMatrix() const
    {
        return invertedLegendreMatrix;
    }


    vector<double> getAlphaVector() const
    {
        return alphaVector;
    }


    std::vector<double> getLagrangeVector() const
    {
        return lagrangeVector;
    }


private: 
    void calculateMesh()
    {
        legendreMesh.resize(2*k);
        std::vector<double> positiveZeros = boost::math::legendre_p_zeros<double>(2*k);

        std::copy(positiveZeros.begin(), positiveZeros.end(), legendreMesh.begin() + k);
        std::transform(positiveZeros.rbegin(), positiveZeros.rend(), legendreMesh.begin(), [](double value){ return -value; });
    }


    void transformMesh()
    {
        transformedMesh.resize(2*k);
        std::transform(legendreMesh.begin(), legendreMesh.end(), transformedMesh.begin(), [this](double value)
        { 
            return transformNode(value); 
        });
    }


    void calculateWeights()
    {
        quadratureWeights.resize(2 * k);

        for (size_t i = 0; i < 2*k; ++i)
        {
            double x_i = legendreMesh[i];

            double P_n = boost::math::legendre_p(2*k, x_i);

            double P_nm1 = boost::math::legendre_p(2*k-1, x_i);

            double dP_n = (2*k)*(P_nm1 - x_i*P_n)/(1.0 - x_i*x_i);
            
            double weight = 2.0 / ((1.0 - x_i * x_i) * dP_n * dP_n);

            double scaled_weight = weight * ((upperBound - lowerBound) / 2.0);

            quadratureWeights[i] = scaled_weight;
        }
    }


    void calculateSquaredAlphas()
    {
        for (size_t i = k; i <= 2*k - 1; ++i)
        {
            measure += (alphaVector[i] * alphaVector[i]);
        }
    }


    void calculateAlphaCoefficients()
    {
        vector<double> lagrangeVectorUblas;
        lagrangeVectorUblas.resize(lagrangeVector.size());

        for (size_t i = 0; i < lagrangeVector.size(); ++i)
        {
            lagrangeVectorUblas(i) = lagrangeVector[i];
        }
        
        alphaVector = prod(invertedLegendreMatrix, lagrangeVectorUblas);
    }


    void evaluateFunctionOnTransformedMesh()
    {
        lagrangeVector.resize(legendreMesh.size());
        std::transform(transformedMesh.begin(), transformedMesh.end(), lagrangeVector.begin(), [this](double value){
            return this->handler.callFunction(value);
        });
    }


    void invertMatrix() 
    {
        matrix<double> copyOfLegendreMatrix(legendreMatrix);
 
        permutation_matrix<std::size_t> pm(copyOfLegendreMatrix.size1());

        int res = lu_factorize(copyOfLegendreMatrix, pm);
        if (res != 0) {
            throw std::invalid_argument("Matrix is singular");
        }

        invertedLegendreMatrix.resize(2*k, 2*k);
        invertedLegendreMatrix.assign(identity_matrix<double>(copyOfLegendreMatrix.size1()));

        lu_substitute(copyOfLegendreMatrix, pm, invertedLegendreMatrix);
    }


    void calculateLegendrePolynomials(matrix<double>& inputMatrix)
    {
        inputMatrix.resize(2*k, 2*k);
        for (size_t i = 0; i < 2*k; ++i)
        {
            for (size_t j = 0; j < legendreMesh.size(); ++j)
            {  
                if (j == 0)
                {
                    inputMatrix(i, j) = 1;
                }
                else
                {
                    inputMatrix(i, j) = transformNode((boost::math::legendre_p(j, legendreMesh[i])));
                }
            }
        }
    }


    static int validateK(int inputK)
    {
        if (inputK > 0)
        {
            return inputK;
        }
        else
        {
            throw std::invalid_argument("k has to be positive");
        }
    }


protected: 
    inline double transformNode(double value) const
    {
        return ((this->upperBound - this->lowerBound) * value + (this->upperBound + this->lowerBound)) / 2;
    }
};

#endif