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
    QuadratureRuleFixture(double lowerBoundInput, double upperBoundInput, const T& handler, double discretizerPrecisionInput = 1e-6, 
    double quadraturePrecisionInput = 1e-4)
    : QuadratureRule<T>(lowerBoundInput, upperBoundInput, handler, discretizerPrecisionInput, quadraturePrecisionInput) {};

    ~QuadratureRuleFixture() {};

    void testCalculateConsolidatedEndpoints()
    {
        this->calculateConsolidatedEndpoints();
    }

    void testEvaluateBasisIntegrals() 
    {
        this->evaluateBasisIntegrals();
    }
};


BOOST_AUTO_TEST_SUITE( QuadratureRuleTestSuite )


BOOST_AUTO_TEST_CASE( TestConstructorPrecisionValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    FunctionHandler<> handlerFunction(testFunction, true);

    double discretizerPrecision = 1e-6;
    double quadraturePrecision = 1e-5;

    BOOST_CHECK_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerFunction, discretizerPrecision, quadraturePrecision), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE( TestConstructorFunctionValidation ) {
    double lowerBound = 10.5;
    double upperBound = 20.4;

    TestClass testObject;
    TestClass* testObjectPtr = &testObject; 
    auto testMethod = std::bind(&TestClass::testMethod, testObjectPtr, std::placeholders::_1);
    std::function<double(const double&)> testMethodPtr = testMethod;

    FunctionHandler<> handlerFunction(testFunction, true);

    FunctionHandler<> handlerMethod(testMethodPtr, true);

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerFunction));

    BOOST_CHECK_NO_THROW(QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerMethod));
}


BOOST_AUTO_TEST_CASE( TestCalculateConsolidatedEndpoints ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5};
    std::vector<double> param2 = {6};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);

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
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
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
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();

    std::vector<double> consolidatedEndpoints = quadrature.getConsolidatedEndpoints();
    std::vector<double> expectedConsolidatedEndpoints = {0, 1, 2};

    BOOST_CHECK_EQUAL(consolidatedEndpoints.size(), 3);
    BOOST_CHECK_EQUAL_COLLECTIONS(consolidatedEndpoints.begin(), consolidatedEndpoints.end(), expectedConsolidatedEndpoints.begin(), expectedConsolidatedEndpoints.end());
}


BOOST_AUTO_TEST_CASE( TestCalculateFinalValues ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();

    std::vector<std::vector<double>> values = quadrature.getValues();

    BOOST_CHECK_EQUAL(values.size(), 4);
    BOOST_CHECK_EQUAL(values[0].size(), 60);
}


BOOST_AUTO_TEST_CASE( TestCalculateFinalWeights) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();

    std::vector<double> weights = quadrature.getWeights();

    BOOST_CHECK_EQUAL(weights.size(), 60);
}


BOOST_AUTO_TEST_CASE( TestCalculateFinalNodes ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();

    std::vector<double> nodes = quadrature.getNodes();

    BOOST_CHECK_EQUAL(nodes.size(), 60);
}


BOOST_AUTO_TEST_CASE( TestCompressorCalling ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();

    BOOST_CHECK_EQUAL(quadrature.getCompressedBasis().size(), 3);
}


BOOST_AUTO_TEST_CASE( TestIntegralEvaluationWithChebyshevNodes ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4, 5, 2, 1};
    std::vector<double> param2 = {6, 3, 9, 8, 7};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    
    std::vector<double> integrals;
    integrals = quadrature.evaluateIntegralsChebyshevNodes();

    displayVector(integrals);
    //TODO calculate exact integarls and have a check. 
}   


BOOST_AUTO_TEST_CASE( TestCompressedBasisSplitting ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();

    std::vector<std::vector<std::vector<double>>> splitCompressedBasis = quadrature.getSplitCompressedBasis();

    BOOST_CHECK_EQUAL(splitCompressedBasis.size(), 3);
    BOOST_CHECK_EQUAL(splitCompressedBasis[0].size(), 2);
    BOOST_CHECK_EQUAL(splitCompressedBasis[0][0].size(), 30);
}


BOOST_AUTO_TEST_CASE( TestCompressedBasisInterpolationCoefficients ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();

    std::vector<std::vector<std::vector<double>>> basisCoefficients = quadrature.getBasisCoefficients();

    BOOST_CHECK_EQUAL(basisCoefficients.size(), 3);
    BOOST_CHECK_EQUAL(basisCoefficients[0].size(), 2);
    BOOST_CHECK_EQUAL(basisCoefficients[0][0].size(), 30);
}


BOOST_AUTO_TEST_CASE( TestBasisIntegralCalculations ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);
    
    QuadratureRuleFixture quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();
    quadrature.testEvaluateBasisIntegrals();

    std::vector<double> basisIntegrals = quadrature.getBasisIntegrals();

    std::cout << "Values of basis integrals " << std::endl;
    displayVector(basisIntegrals);

    BOOST_CHECK_EQUAL(basisIntegrals.size(), 3);
}


BOOST_AUTO_TEST_SUITE_END()