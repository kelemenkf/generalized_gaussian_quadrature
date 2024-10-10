#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include "ggq.hpp"
using namespace boost::numeric::ublas;


template<typename InputClass>
class Discretizer: public QuadratureRule<InputClass>
{
private: 
    int k;
    double precision;
    std::vector<double> legendreMesh;
    std::vector<double> transformedMesh;
    matrix<double> legendreMatrix;
    matrix<double> invertedLegendreMatrix;
    std::vector<double> lagrangeVector;
    std::vector<double> alphaVector;

    using InputMethodType = double (InputClass::*)(const double&);
    using InputFunctionType = double(*)(const double&);


public:
    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : QuadratureRule<InputClass>(lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr), 
    k(validateK(kInput)), precision(validatePrecision(precisionInput)) {
        calculateMesh();
        transformMesh();
        calculateLegendrePolynomials(legendreMatrix);
        invertMatrix();
    };

    ~Discretizer() {};

    std::vector<double> discretizationRoutine()
    {
        evaluateFunctionOnTransformedMesh();
        //After a single funciton is discretized the object should clear these, so this function can be called again, or recursivley
        // lagrangeVector.clear();
        // alphaVector.clear();
    }

private: 
    void evaluateFunctionOnTransformedMesh()
    {
        lagrangeVector.resize(legendreMesh.size());
        std::transform(lagrangeVector.begin(), lagrangeVector.end(), lagrangeVector.begin(), [this](double value){
            return this->functionPtr(value);
        });
    }


    void invertMatrix() {
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
                if (i == 0)
                {
                    inputMatrix(i, j) = 1;
                }
                else
                {
                    inputMatrix(i, j) = transformNode((boost::math::legendre_p(i, legendreMesh[j])));
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

    static int validatePrecision(double inputPrecision)
    {
        if (inputPrecision > 0)
        {
            return inputPrecision;
        }
        else 
        {
            throw std::invalid_argument("Precision has to be positive");
        }
    }


protected: 
    inline double transformNode(double value) const
    {
        return ((this->upperBound - this->lowerBound) * value + (this->upperBound + this->lowerBound)) / 2;
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
};

#endif