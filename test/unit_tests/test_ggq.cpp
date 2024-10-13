#define BOOST_TEST_MODULE QuadratureRuleTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "ggq.hpp"


class TestClass
{
public:
    double testMethod(const double& x)
    {
        return x;
    }
};


double testFunction(const double& x)
{
    return x;
}


struct QuadratureRuleFixture: public QuadratureRule
{
    using InputFunctionType = std::function<double(const double&)>;

    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, InputFunctionType function) : QuadratureRule(lowerBoundInput, upperBoundInput, function) {};

    ~QuadratureRuleFixture() {};
};


BOOST_AUTO_TEST_SUITE( QuadratureRuleTestSuite )


BOOST_AUTO_TEST_CASE( TestConstructorValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    TestClass testObject;
    TestClass* testObjectPtr = &testObject; 
    auto testMethod = std::bind(&TestClass::testMethod, &testObject, std::placeholders::_1);

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, testFunction));

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, testMethod));
}


BOOST_AUTO_TEST_SUITE_END()