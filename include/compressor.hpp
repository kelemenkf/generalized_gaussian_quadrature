#ifndef COMPRESSOR_HPP
#define COMPRESSOR_HPP

#include "ggq.hpp"


template<typename T>
class Compressor
{
private:
    T quadrature;
    matrix<double> A;
    matrix<double> U;
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


    matrix<double> getA() const
    {
        return A;
    }


private:
    void constructA()
    {
        A.resize(nodes.size(), values.size());
        for (size_t column = 0; column < values.size(); ++column)
        {
            for (size_t row = 0; row < nodes.size(); ++row)
            {
                A(row, column) = values[column][row];
            }
        }
    }


    void decomposeIntoQR()
    {

    }


    void scaleU()
    {

    }

protected:
};

#endif