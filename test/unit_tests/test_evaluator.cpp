#define BOOST_TEST_MODULE EvaluatorTestSuite 
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "evaluator.hpp"
#include "ggq.hpp"
#include "utils.hpp"


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
FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);

QuadratureRule<FunctionHandler<std::vector<double>, std::vector<double>>> quadratureObject(lowerBound, upperBound, handlerPiecewiseSmooth);


BOOST_AUTO_TEST_SUITE( EvaluatorTestSuite )


BOOST_AUTO_TEST_CASE( TestEvaluatorConstructor ) {
    quadratureObject.calculateQuadratureNodes();
    quadratureObject.compressFunctionSpace();
    std::vector<std::vector<double>> basisFunctions = quadratureObject.getCompressedBasis();
    std::vector<double> nodes = quadratureObject.getNodes();
    std::vector<double> endpoints = quadratureObject.getConsolidatedEndpoints();

    Evaluator evaluator(nodes, basisFunctions[0], endpoints);

    std::vector<vector<double>> coefficients = evaluator.getCoefficients();

    size_t intervalLength = 30;

    // for (size_t i = 0; i < coefficients.size(); ++i)
    // {
    //     for (size_t j = 0; j < intervalLength; ++j)
    //     {
    //         std::cout << i << " " << intervalLength * i + j << std::endl;
    //         BOOST_CHECK_CLOSE_FRACTION(coefficients[i].evaluate(nodes[intervalLength * i + j]), basisFunctions[0][intervalLength * i + j], 1e-9);
    //     }     
    // }

    BOOST_CHECK_CLOSE_FRACTION(evaluator.evaluate(nodes[0], coefficients[0]), basisFunctions[0][0], 1e-9);
}


BOOST_AUTO_TEST_CASE( TestMatrixInversion ) {
    
}

BOOST_AUTO_TEST_SUITE_END()