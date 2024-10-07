#include "ggq.hpp"

template<typename InputClass>
QuadratureRule<InputClass>::QuadratureRule(InputFunctionType function, InputClass& inputClass, double lowerBoundInput, double upperBoundInput) : 
lowerBound(lowerBoundInput), upperBound(upperBoundInput), functionPtr(function), objectRef(inputClass)
{

};


template<typename InputClass>
QuadratureRule<InputClass>::~QuadratureRule() 
{

};


template<typename InputClass>
double QuadratureRule<InputClass>::validateLowerBound(double input)
{
    //TODO does this have a bound?
}


template<typename InputClass>
double QuadratureRule<InputClass>::validateUpperBound(double input)
{
    //TODO does this have a bound?
}


