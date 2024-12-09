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

    displayVector(basisFunctions[0]);
    displayVector(nodes);

    Evaluator evaluator(nodes, basisFunctions[0]);

    polynomial<double> coefficients = evaluator.getCoefficients();

    std::cout << coefficients << std::endl;

    std::cout << coefficients.evaluate(nodes[0]) << std::endl;

    BOOST_CHECK_EQUAL(coefficients.evaluate(nodes[0]), basisFunctions[0][0]);
}





BOOST_AUTO_TEST_CASE( TestEvaluatorSimplePolynomial ) {
    quadratureObject.calculateQuadratureNodes();
    quadratureObject.compressFunctionSpace();
    std::vector<double> values{0.3, 0.5, 0.6, 1.2};
    std::vector<double> legendreMesh(4);
    std::vector<double> positiveZeros = boost::math::legendre_p_zeros<double>(4);
    std::copy(positiveZeros.begin(), positiveZeros.end(), legendreMesh.begin() + 2);
    std::transform(positiveZeros.rbegin(), positiveZeros.rend(), legendreMesh.begin(), [](double value){ return -value; });

    Evaluator evaluator(legendreMesh, values);

    displayVector(positiveZeros);

    polynomial<double> coefficients = evaluator.getCoefficients();

    std::cout << coefficients << std::endl;

    std::cout << coefficients.evaluate(legendreMesh[0]) << std::endl;

    BOOST_CHECK_CLOSE(coefficients.evaluate(legendreMesh[0]), values[0], 1e-9);
}



BOOST_AUTO_TEST_SUITE_END()