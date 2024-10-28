#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP


#include <boost/math/special_functions/legendre.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
using namespace boost::numeric::ublas;


template<typename T>
class Interpolator 
{
private:
    int k;
    double lowerBound;
    double upperBound;
    std::vector<double> legendreMesh;
    std::vector<double> transformedMesh;
    matrix<double> legendreMatrix;
    matrix<double> invertedLegendreMatrix;
    vector<double> alphaVector;
    std::vector<double> lagrangeVector; 
    T handler;


public: 
    Interpolator(int kInput, double lowerBoundInput, double upperBoundInput, std::vector<double> transformedMeshInput, std::vector<double> legendreMeshInput) : 
    k(kInput), lowerBound(lowerBoundInput), upperBound(upperBoundInput), transformedMesh(transformedMeshInput), legendreMesh(legendreMeshInput)
    {

    }


    ~Interpolator() 
    {

    }


    void interpolatePolynomial()
    {
        calculateLegendrePolynomials(legendreMatrix);
        invertMatrix();
        evaluateFunctionOnTransformedMesh();
        calculateAlphaCoefficients();
    }

private:
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
        lagrangeVector.resize(transformedMesh.size());
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
            for (size_t j = 0; j < transformedMesh.size(); ++j)
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


protected:
    inline double transformNode(double value) const
    {
        return ((this->upperBound - this->lowerBound) * value + (this->upperBound + this->lowerBound)) / 2;
    }
};


#endif