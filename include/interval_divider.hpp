#ifndef INTERVAL_DIVIDER_HPP
#define INTERVAL_DIVIDER_HPP

#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <cmath>
using namespace boost::numeric::ublas;


template<typename T>
class IntervalDivider
{
private: 
    int k;
    std::vector<double> legendreMesh;
    std::vector<double> transformedMesh;
    matrix<double> legendreMatrix;
    matrix<double> invertedLegendreMatrix;
    vector<double> alphaVector;
    std::vector<double> lagrangeVector;
    double measure;
    double lowerBound;
    double upperBound;
    T handler;

public: 
    IntervalDivider(int kInput, double lowerBoundInput, double upperBoundInput, const T& handlerInput) 
    : k(kInput), lowerBound(lowerBoundInput), upperBound(upperBoundInput), handler(handlerInput)
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
        evaluateFunctionOnTransformedMesh();
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


    matrix<double> getLegendreMatrix() const
    {
        return legendreMatrix;
    }


    matrix<double> getInvertedLegendreMatrix() const
    {
        return invertedLegendreMatrix;
    }


    std::vector<double> getLagrangeVector() const
    {
        return lagrangeVector;
    }


    vector<double> getAlphaVector() const
    {
        return alphaVector;
    }


private: 
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