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
    MatrixXd A;
    MatrixXd U;
    MatrixXd scaledU;
    MatrixXd R;
    std::vector<double> nodes;
    std::vector<double> weights;
    std::vector<std::vector<double>> values;
    std::vector<double> normalizingFactors;


public:
    Compressor(const T& quadratureInput) : quadrature(quadratureInput)
    {
        nodes = quadrature.getNodes();
        weights = quadrature.getWeights();
        values = quadrature.getValues();
        constructA();
        decomposeIntoQR();
        scaleU();
        getNormalizingFactors();
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
        std::cout << "Shape of R " << R.rows() << "x" << R.cols() << std::endl;
    }


    void scaleU()
    {
        scaledU = U;
        for (size_t row = 0; row < weights.size(); ++row)
        {
            scaledU.row(row) /= sqrt(weights[row]);
        }
    }


    void getNormalizingFactors()
    {  
        for (size_t i = 0; i < R.rows(); ++i)
        {
            for (size_t j = 0; j < R.cols(); ++j)
            {
                if (i == j)
                {
                    std::cout << R(i,j) << std::endl;
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