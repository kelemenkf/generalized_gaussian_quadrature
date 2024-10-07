#define BOOST_TEST_MODULE QuadratureRuleTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "ggq.hpp"


class TestClass
{
public:
    double testMethod(const std::vector<double>& input)
    {
        std::cout << "Test method ran" << std::endl;
        return input[0];
    }
};


double testFunction(const std::vector<double>& input)
{
    std::cout << "Test function ran" << std::endl;
    return input[0];
}


template<typename InputClass = double>
struct QuadratureRuleFixture: public QuadratureRule<InputClass>
{
    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, InputTypeFunction function = NULL, 
    InputTypeMethod method = NULL, InputClass& inputObject = NULL) : QuadratureRule(lowerBoundInput, upperBoundInput, function, method, 
    inputObject) {};

    ~QuadratureRuleFixture() {};

    InputFunctionType* testGetFunctionPointer()
    {
        return getFunctionPointer();
    }

    InputMethodType* testGetMethodPointer()
    {
        return getMethodPointer();
    }
};


BOOST_AUTO_TEST_SUITE( QuadratureRuleTestSuite )


BOOST_AUTO_TEST_CASE( TestConstructorValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    QuadratureRuleFixture quadrature(lowerBound, upperBound);
}


BOOST_AUTO_TEST_SUITE_END()