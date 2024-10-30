#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include <Eigen/Dense>
#include <numeric>
#include "utils.hpp"
using namespace Eigen;


template<typename T>
class Compressor
{
private:
    const T* quadrature;
    double quadraturePrecision;
    MatrixXd A;
    MatrixXd U;
    MatrixXd scaledU;
    MatrixXd R;
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;
    std::vector<double> normalizingFactors;
    std::vector<std::vector<double>> scaledDiscardedU;
    std::vector<double> rVector;
    MatrixXd B;


public:
    Compressor(const T* quadratureInput, double quadraturePrecisionInput = 1e-3) : quadrature(quadratureInput), quadraturePrecision(quadraturePrecisionInput)
    {
        nodes = quadrature->getNodes();
        weights = quadrature->getWeights();
        values = quadrature->getValues();
        constructA();
        decomposeIntoQR();
        scaleU();
        calculateNormalizingFactors();
        discardFunctions(); 
        calculateRVector();
        constructB();
    }


    ~Compressor() 
    {

    }


    MatrixXd getA() const
    {
        return A;
    }


    MatrixXd getU() const
    {
        return U;
    }


    MatrixXd getScaledU() const
    {
        return scaledU;
    }


    std::vector<double> getNormalizingFactors() const 
    {
        return normalizingFactors;
    }


    std::vector<std::vector<double>> getCompressedBasis() const
    {
        return scaledDiscardedU;
    }


    std::vector<double> getRVector() const
    {
        return rVector;
    }


    MatrixXd getB() const
    {
        return B;
    }


private:
    void constructA()
    {
        A.resize(nodes.size(), values.size());
        for (size_t column = 0; column < values.size(); ++column)
        {
            for (size_t row = 0; row < nodes.size(); ++row)
            {
                A(row, column) = values[column][row] * sqrt(weights[row]);
            }
        }
    }


    void decomposeIntoQR()
    {
        HouseholderQR<MatrixXd> qr(A);

        U = qr.householderQ();
        R = qr.matrixQR().triangularView<Eigen::Upper>();
    }


    void scaleU()
    {
        scaledU = U;
        for (size_t row = 0; row < weights.size(); ++row)
        {
            scaledU.row(row) /= sqrt(weights[row]);
        }
    }


    void calculateNormalizingFactors()
    {  
        for (size_t i = 0; i < R.rows(); ++i)
        {
            for (size_t j = 0; j < R.cols(); ++j)
            {
                if (i == j)
                {
                    if (R(i,j) < 0)
                    {
                        //Correction for Householder QR producing negative factors
                        normalizingFactors.push_back(-R(i,j));
                    }
                    else 
                    {
                        normalizingFactors.push_back(R(i,j));
                    }
                }
            }
        }
    }


    void discardFunctions()
    {
        for (size_t i = 0; i < normalizingFactors.size(); ++i)
        {
            if (normalizingFactors[i] > quadraturePrecision)
            {
                Eigen::VectorXd column = scaledU.col(i);
                std::vector<double> vec(column.data(), column.data() + column.size());
                scaledDiscardedU.push_back(vec);
            }
        }
    }


    double innerProduct(size_t index)
    {
        double result = 0; 

        if (scaledDiscardedU[index].size() != weights.size()) {
            throw std::runtime_error("Size mismatch between vectors.");
        }

        for (size_t i = 0; i < scaledDiscardedU[index].size(); ++i)
        {
            result += (scaledDiscardedU[index][i] * weights[i]);
        }

        return result;
    }


    void calculateRVector()
    {
        for (size_t i = 0; i < scaledDiscardedU.size(); ++i)
        {
            double res = innerProduct(i);
            rVector.push_back(res);
        }
    }


    void constructB()
    {
        B.resize(scaledDiscardedU.size(), weights.size());
        for (size_t row = 0; row < scaledDiscardedU.size(); ++row)
        {
            for (size_t column = 0; column < weights.size(); ++column)
            {
                B(row, column) = scaledDiscardedU[row][column] * sqrt(weights[column]);
            }
        }
    }
};
 

#endif