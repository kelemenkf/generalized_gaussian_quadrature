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

template<typename T>
struct QuadratureRuleFixture: public QuadratureRule<T>
{
    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, const T& handler)
    : QuadratureRule<T>(lowerBoundInput, upperBoundInput, handler) {};

    ~QuadratureRuleFixture() {};
};


BOOST_AUTO_TEST_SUITE( QuadratureRuleTestSuite )


BOOST_AUTO_TEST_CASE( TestConstructorValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    TestClass testObject;
    TestClass* testObjectPtr = &testObject; 
    auto testMethod = std::bind(&TestClass::testMethod, testObjectPtr, std::placeholders::_1);
    std::function<double(const double&)> testMethodPtr = testMethod;

    FunctionHandler<> handlerFunction(testFunction);

    FunctionHandler<> handlerMethod(testMethodPtr);

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerFunction));

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerMethod));
}


BOOST_AUTO_TEST_SUITE_END()