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


template<typename InputClass>
struct QuadratureRuleFixture: public QuadratureRule<InputClass>
{
    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    using InputFunctionType = double (*)(const std::vector<double>&);

    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, InputFunctionType function = nullptr, 
    InputMethodType method = nullptr, InputClass* inputObject = nullptr) : QuadratureRule<InputClass>(lowerBoundInput, upperBoundInput, function, method, 
    inputObject) {};

    ~QuadratureRuleFixture() {};
};


BOOST_AUTO_TEST_SUITE( QuadratureRuleTestSuite )


BOOST_AUTO_TEST_CASE( TestConstructorValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    TestClass testObject;

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture<TestClass> quadrature(lowerBound, upperBound, &testFunction));

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture<TestClass> quadrature(lowerBound, upperBound, nullptr, &TestClass::testMethod, &testObject));

    BOOST_CHECK_THROW(QuadratureRuleFixture<TestClass> quadrature(lowerBound, upperBound), std::invalid_argument);
}


BOOST_AUTO_TEST_SUITE_END()