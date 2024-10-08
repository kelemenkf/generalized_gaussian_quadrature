#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include <boost/math/special_functions/legendre.hpp>
#include "ggq.hpp"


template<typename InputClass>
class Discretizer: public QuadratureRule<InputClass>
{
private: 
    int k;
    double precision;
    std::vector<double> legendreMesh;

    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    using InputFunctionType = double(*)(const std::vector<double>&);


public:
    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : QuadratureRule<InputClass>(lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr), 
    k(validateK(kInput)), precision(validatePrecision(precisionInput)) {
        calculateMesh();

        for (auto element = legendreMesh.begin(); element != legendreMesh.end(); ++element)
        {
            std::cout << *element << " ";
        }
        std::cout << std::endl;
    };

    ~Discretizer() {};


private: 
    void calculateMesh()
    {
        legendreMesh.resize(2*k);
        std::vector<double> positiveZeros = boost::math::legendre_p_zeros<double>(2*k);

        std::copy(positiveZeros.begin(), positiveZeros.end(), legendreMesh.begin() + k);
        std::transform(positiveZeros.rbegin(), positiveZeros.rend(), legendreMesh.begin(), [](double value){ return -value; });
    }

    void transformMesh();

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
    std::vector<double> getLegendreMesh() const
    {
        return legendreMesh;
    }
};

#endif