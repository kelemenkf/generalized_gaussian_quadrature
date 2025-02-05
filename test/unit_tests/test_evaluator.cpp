#define BOOST_TEST_MODULE EvaluatorTestSuite 
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "evaluator.hpp"
#include "ggq.hpp"


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



double lowerBound = 0;
double upperBound = 2; 
std::vector<double> param1 = {5, 4};
std::vector<double> param2 = {6, 3};
FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);

QuadratureRule<FunctionHandler<std::vector<double>, std::vector<double>>> quadratureObject(lowerBound, upperBound, handlerPiecewiseSmooth);


BOOST_AUTO_TEST_SUITE( EvaluatorTestSuite )


BOOST_AUTO_TEST_CASE( TestEvaluatorConstructor ) {
    quadratureObject.calculateQuadratureNodes();
    quadratureObject.compressFunctionSpace();
    std::vector<std::vector<double>> basisFunctions = quadratureObject.getCompressedBasis();
    std::vector<double> nodes = quadratureObject.getNodes();
    std::vector<double> endpoints = quadratureObject.getConsolidatedEndpoints();
}


BOOST_AUTO_TEST_CASE( TestNodeInversion ) {
    std::vector<std::vector<double>> basisFunctions = quadratureObject.getCompressedBasis();
    std::vector<double> endpoints = quadratureObject.getConsolidatedEndpoints();
    std::vector<double> y = quadratureObject.getBasisFunctionInterval(0, 0);
    std::vector<std::vector<double>> splitNodes = quadratureObject.getSplitNodes();

    Interpolator divider(15, endpoints[0], endpoints[1], handlerPiecewiseSmooth, y);
    divider.interpolateFunction();

    std::vector<double> alphas = divider.getAlphaVector();

    Evaluator evaluator(alphas, endpoints[0], endpoints[1], splitNodes[0]);
    evaluator.evaluateInput();

    std::vector<double> reversedNodes = evaluator.getReversedNodes();
    std::vector<double> legendreNodes = divider.getLegendreMesh();

    for (size_t i = 0; i < reversedNodes.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(reversedNodes[i], legendreNodes[i], 1e-9);
    }
}


BOOST_AUTO_TEST_CASE( TestMatrixInversion ) {
    std::vector<std::vector<double>> basisFunctions = quadratureObject.getCompressedBasis();
    std::vector<double> nodes = quadratureObject.getNodes();
    std::vector<double> endpoints = quadratureObject.getConsolidatedEndpoints();
    std::vector<double> y = quadratureObject.getBasisFunctionInterval(0, 0);

    Interpolator divider(15, endpoints[0], endpoints[1], handlerPiecewiseSmooth, y);

    divider.interpolateFunction();

    std::vector<double> alphaVector = divider.getAlphaVector(); 
    matrix<double> legendreMatrix = divider.getLegendreMatrix();
    matrix<double> invertedLegendreMatrix = divider.getInvertedLegendreMatrix();
    std::vector<double> lagrangeVector = divider.getLagrangeVector();

    vector<double> calculatedF = prod(legendreMatrix, convertStdVectorToBoost(alphaVector));  

    for (size_t i = 0; i < lagrangeVector.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(calculatedF(i), lagrangeVector[i], 1e-9);
    }
}


BOOST_AUTO_TEST_CASE( TestNodeEvaluation ) {
    std::vector<std::vector<double>> basisFunctions = quadratureObject.getCompressedBasis();
    std::vector<double> endpoints = quadratureObject.getConsolidatedEndpoints();
    std::vector<double> y = quadratureObject.getBasisFunctionInterval(0, 0);
    std::vector<std::vector<double>> splitNodes = quadratureObject.getSplitNodes();

    Interpolator divider(15, endpoints[0], endpoints[1], handlerPiecewiseSmooth, y);
    divider.interpolateFunction();

    std::vector<double> alphas = divider.getAlphaVector();

    Evaluator evaluator(alphas, endpoints[0], endpoints[1], splitNodes[0]);
    evaluator.evaluateInput();

    std::vector<double> reversedNodes = evaluator.getReversedNodes();
    std::vector<double> legendreNodes = divider.getLegendreMesh();

    for (size_t i = 0; i < reversedNodes.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(reversedNodes[i], legendreNodes[i], 1e-9);
    }
}



BOOST_AUTO_TEST_CASE(TestSmallerNumberOfNodes ) {
}


BOOST_AUTO_TEST_CASE( TestFirstDerivativeEvaluator ) {
    std::vector<double> coefficients{1, 2, 3};
    std::vector<double> inputNodes{1, 2, 3};

    Evaluator evaluator(coefficients, 2, 3, inputNodes);

    double result = evaluator.evaluateFirstDerivative(2); 

    std::cout << result << std::endl;

    BOOST_CHECK_EQUAL(result, 11);
}


BOOST_AUTO_TEST_SUITE_END()