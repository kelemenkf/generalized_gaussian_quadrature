#ifndef DISCRETIZER_HPP
#define DISCRETIZER_HPP

#include "ggq.hpp"


template<typename InputClass>
class Discretizer: public QuadratureRule<InputClass>
{
private: 
    int k;
    double precision;

    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    using InputFunctionType = double(*)(const std::vector<double>&);


public:
    Discretizer(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : QuadratureRule<InputClass>(lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr), 
    k(validateK(kInput)), precision(validatePrecision(precisionInput)) {};

    ~Discretizer() {};


private: 
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
};

#endif