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

    }


    void decomposeIntoQR()
    {

    }


    void scaleU()
    {

    }

protected:
};