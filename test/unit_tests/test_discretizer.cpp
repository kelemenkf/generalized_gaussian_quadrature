#define BOOST_TEST_MODULE DiscretizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "discretizer.hpp"


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


template<typename InputClass>
struct DiscretizerFixture: public Discretizer<InputClass>
{
    using InputMethodType = double (InputClass::*)(const double&);
    using InputFunctionType = double(*)(const double&);

    DiscretizerFixture(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : Discretizer<InputClass>(kInput, precisionInput, lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr) {};

    ~DiscretizerFixture() {};

    double testGetPrecision() const
    {
        return this->getPrecision();
    }

    bool testEvaluateStoppingCondition()
    {
        return this->evaluateStoppingCondition();
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

    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(-10, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(30, -0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_NO_THROW(Discretizer<TestClass> discretizer(30, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr));
}


// BOOST_AUTO_TEST_CASE( TestEvaluateStoppingCondition ) {
//     double lowerBound = 1;
//     double upperBound = 2;
//     const int k = 30;
//     const double precision = 1e-6;

//     DiscretizerFixture<TestClass> discretizer(k, precision, lowerBound, upperBound, singularTestFunction);
//     discretizer.discretizationRoutine();
// }


BOOST_AUTO_TEST_CASE( TestDiscretizerCalculateMeasures ) {
    double lowerBound = 1;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture<TestClass> discretizer(k, precision, lowerBound, upperBound, testFunction);
    discretizer.testCalculateMeasures();

    std::vector<double> measures = discretizer.testGetMeasureVector();

    BOOST_CHECK_EQUAL(measures.size(), 1);
}
 

BOOST_AUTO_TEST_CASE( TestDiscretizerFindEndpoints ) {
    double lowerBound = 1;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture<TestClass> discretizer(k, precision, lowerBound, upperBound, testFunction);
    discretizer.discretizationRoutine();

    std::vector<double> endpoints = discretizer.getFinalEndpoints();
    std::vector<double> expectedEndpoints = {lowerBound, upperBound};

    BOOST_CHECK_EQUAL(endpoints.size(), expectedEndpoints.size());

    for (size_t i = 0; i < endpoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(endpoints[i], expectedEndpoints[i], 1e-6);
    }
}


BOOST_AUTO_TEST_CASE( TestDiscretizerFindEndpointsSingularFunction ) {
    double lowerBound = 1;
    double upperBound = 2;
    const int k = 30;
    const double precision = 1e-6;

    DiscretizerFixture<TestClass> discretizer(k, precision, lowerBound, upperBound, singularTestFunction);
    discretizer.discretizationRoutine();

    std::vector<double> endpoints = discretizer.getFinalEndpoints();
    std::vector<double> expectedEndpoints = {lowerBound, 0, upperBound};

    BOOST_CHECK_EQUAL(endpoints.size(), expectedEndpoints.size());

    for (size_t i = 0; i < endpoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(endpoints[i], expectedEndpoints[i], 1e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END()