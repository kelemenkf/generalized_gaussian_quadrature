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


public:
    Compressor(const T& quadratureInput) : quadrature(quadratureInput)
    {
        nodes = quadrature.getNodes();
        weights = quadrature.getWeights();
        values = quadrature.getValues();
        constructA();
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
                A(row, column) = values[column][row] * sqrt(weights[column]);
            }
        }
    }


    void decomposeIntoQR()
    {
        HouseholderQR<MatrixXd> qr(A);

        U = qr.householderQ();
        R = qr.matrixQR().triangularView<Eigen::Upper>();

        std::cout << U << std::endl;
        std::cout << R << std::endl;
    }


    void scaleU()
    {
        scaledU = U;
        for (size_t column = 0; column < values.size(); ++column)
        {
            scaledU.col(column) /= sqrt(weights[column]);
        }
    }


protected:
};

#endif