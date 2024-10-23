#include "function_handler.hpp"
#include "interval_divider.hpp"

template<typename T>
class Compressor
{
private:
    T handler;
    std::vector<double> weights;
    std::vector<double> nodes;
    std::vector<std::vector<double>> values;
    matrix<double> A;
    matrix<double> U;


public:
    Compressor() 
    {

    }

    ~Compressor() 
    {

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