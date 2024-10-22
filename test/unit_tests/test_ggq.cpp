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


double testFunction2ParamPC(const double& x, const double& param1, const double& param2)
{
    if (x <= 1)
    {
        return param1 * x * x;
    }
    else 
    {
        return param2 - x;
    }
}


template<typename T>
struct QuadratureRuleFixture: public QuadratureRule<T>
{
    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, const T& handler)
    : QuadratureRule<T>(lowerBoundInput, upperBoundInput, handler) {};

    ~QuadratureRuleFixture() {};

    void testCalculateConsolidatedEndpoints()
    {
        this->calculateConsolidatedEndpoints();
    }
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


BOOST_AUTO_TEST_CASE( TestCalculateConsolidatedEndpoints ) {
   double lowerBound = 0;
   double upperBound = 2; 
   int k = 30;
    std::vector<double> param1 = {5};
    std::vector<double> param2 = {6};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);

   QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
   quadrature.testCalculateConsolidatedEndpoints();

   std::vector<double> consolidatedEndpoints = quadrature.getConsolidatedEndpoints();
   std::vector<double> expectedConsolidatedEndpoints = {0, 1, 2};

    BOOST_CHECK_EQUAL(consolidatedEndpoints.size(), 3);
    BOOST_CHECK_EQUAL_COLLECTIONS(consolidatedEndpoints.begin(), consolidatedEndpoints.end(), expectedConsolidatedEndpoints.begin(), expectedConsolidatedEndpoints.end());
}


BOOST_AUTO_TEST_CASE( TestCalculateConsolidatedEndpointsMoreParameters ) {
   double lowerBound = 0;
   double upperBound = 2; 
   int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);
    
   QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
   quadrature.testCalculateConsolidatedEndpoints();

   std::vector<double> consolidatedEndpoints = quadrature.getConsolidatedEndpoints();
   std::vector<double> expectedConsolidatedEndpoints = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};

    BOOST_CHECK_EQUAL(consolidatedEndpoints.size(), 12);
    BOOST_CHECK_EQUAL_COLLECTIONS(consolidatedEndpoints.begin(), consolidatedEndpoints.end(), expectedConsolidatedEndpoints.begin(), expectedConsolidatedEndpoints.end());
}



BOOST_AUTO_TEST_CASE( TestCalculateFinalConsolidatedEndpointsMoreParameters ) {
   double lowerBound = 0;
   double upperBound = 2; 
   int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);
    
   QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
   quadrature.calculateQuadratureNodes();

   std::vector<double> consolidatedEndpoints = quadrature.getConsolidatedEndpoints();
   std::vector<double> expectedConsolidatedEndpoints = {0, 1, 2};

    BOOST_CHECK_EQUAL(consolidatedEndpoints.size(), 3);
    BOOST_CHECK_EQUAL_COLLECTIONS(consolidatedEndpoints.begin(), consolidatedEndpoints.end(), expectedConsolidatedEndpoints.begin(), expectedConsolidatedEndpoints.end());
}


BOOST_AUTO_TEST_SUITE_END()