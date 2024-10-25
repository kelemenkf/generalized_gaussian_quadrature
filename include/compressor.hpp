#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include "ggq.hpp"
#include <Eigen/Dense>
using namespace Eigen;


template<typename T>
class Compressor
{
private:
    T quadrature;
    double quadraturePrecision;
    MatrixXd A;
    MatrixXd U;
    MatrixXd scaledU;
    MatrixXd R;
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;
    std::vector<double> normalizingFactors;


public:
    Compressor(const T& quadratureInput, double quadraturePrecisionInput = 1e-3) : quadrature(quadratureInput), quadraturePrecision(quadraturePrecisionInput)
    {
        nodes = quadrature.getNodes();
        weights = quadrature.getWeights();
        values = quadrature.getValues();
        constructA();
        decomposeIntoQR();
        scaleU();
        calculateNormalizingFactors();
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
                    normalizingFactors.push_back(R(i,j));
                }
            }
        }
    }


    void discardFunctions()
    {

    }


protected:
};

#endif