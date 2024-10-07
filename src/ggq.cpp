#include "ggq.hpp"

template<typename InputClass>
QuadratureRule<InputClass>::QuadratureRule(double lowerBoundInput, double upperBoundInput, InputFunctionType function = NULL, 
InputMethodType methd = NULL, InputClass& inputClass = NULL) : 
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


