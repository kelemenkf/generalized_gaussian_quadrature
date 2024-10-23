#define BOOST_TEST_MODULE DiscretizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "discretizer.hpp"
#include "function_handler.hpp"


class TestClass
{
public:
    double testMethod(const double& x)
    {
        return (5*std::pow(x,3) - 3*x) / 2;
    }
};


double testFunction(const double& x)
{
    return (5*std::pow(x,3) - 3*x) / 2;
}


double singularTestFunction(const double& x)
{
    return 1 / x;
}


double piecewiseSmoothFunction(const double& x)
{
    if (x <= 1)
    {
        return x * x;
    }
    else 
    {
        return 2 - x;
    }
}


double highlyOscillatoryFunction(const double& x)
{
    return sin(50 * x);
}


FunctionHandler<> handlerPolynomial(testFunction);
FunctionHandler<> handlerPiecewiseSmooth(piecewiseSmoothFunction);
FunctionHandler<> handlerHighlyOscillatory(highlyOscillatoryFunction);


template<typename T>
struct DiscretizerFixture: public Discretizer<T>
{
    DiscretizerFixture(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, const T& handlerInput) 
    : Discretizer<T> (kInput, precisionInput, lowerBoundInput, upperBoundInput, handlerInput) {}

    ~DiscretizerFixture() {}

    double testGetPrecision() const
    {
        return this->getPrecision();
    }

    void testCalculateMeasures()
    {
        this->calculateMeasures();
    }

    std::vector<double> testGetMeasureVector() const
    {
        return this->getMeasureVector();
    }
};


BOOST_AUTO_TEST_SUITE( DiscretizerTestSuite )

BOOST_AUTO_TEST_CASE( TestDiscretizerValidation ) {
    double lowerBound = 1;
    double upperBound = 2;

    BOOST_CHECK_THROW(DiscretizerFixture discretizer(-10, 0.01, lowerBound, upperBound, handlerPolynomial), std::invalid_argument);
    BOOST_CHECK_THROW(DiscretizerFixture discretizer(30, -0.01, lowerBound, upperBound, handlerPolynomial), std::invalid_argument);
    BOOST_CHECK_NO_THROW(DiscretizerFixture discretizer(30, 0.01, lowerBound, upperBound, handlerPolynomial));
}


BOOST_AUTO_TEST_CASE( TestDiscretizerCalculateMeasures ) {
    double lowerBound = 1;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture discretizer(k, precision, lowerBound, upperBound, handlerPolynomial);
    
    std::vector<double> measures = discretizer.testGetMeasureVector();

    BOOST_CHECK_EQUAL(measures.size(), 1);
}
 

BOOST_AUTO_TEST_CASE( TestDiscretizerFindEndpoints ) {
    double lowerBound = 1;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture discretizer(k, precision, lowerBound, upperBound, handlerPolynomial);

    std::vector<double> endpoints = discretizer.getFinalEndpoints();
    std::vector<double> expectedEndpoints = {lowerBound, upperBound};

    BOOST_CHECK_EQUAL(endpoints.size(), expectedEndpoints.size());

    for (size_t i = 0; i < endpoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(endpoints[i], expectedEndpoints[i], 1e-6);
    }
}


BOOST_AUTO_TEST_CASE( TestDiscretizerFindEndpointsPiecewiseSmoothFunction ) {
    double lowerBound = 0;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture discretizer(k, precision, lowerBound, upperBound, handlerPiecewiseSmooth);

    std::vector<double> endpoints = discretizer.getFinalEndpoints();
    std::vector<double> expectedEndpoints = {lowerBound, 1, upperBound};

    BOOST_CHECK_EQUAL(endpoints.size(), expectedEndpoints.size());

    for (size_t i = 0; i < endpoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(endpoints[i], expectedEndpoints[i], 1e-6);
    }
}


BOOST_AUTO_TEST_SUITE_END()