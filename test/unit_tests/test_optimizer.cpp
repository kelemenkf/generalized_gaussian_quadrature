#define BOOST_TEST_MODULE OptimizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
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

struct OptimizerFixture : public Optimizer
{
    OptimizerFixture(const std::vector<double>& inputChebyshevNodes, const std::vector<double>& inputChebyshevWeights, 
    const std::vector<std::vector<std::vector<double>>>& inputBasisCoefficients, const std::vector<std::vector<std::vector<double>>>& 
    inputSplitCompressedBasis, const std::vector<double>& inputBasisIntegrals)
    : Optimizer(inputChebyshevNodes, inputChebyshevWeights, inputBasisCoefficients, inputSplitCompressedBasis, inputBasisIntegrals) 
    {};

    ~OptimizerFixture() {};

    void testFormJacobian() 
    {
        this->formJacobian();
    }
};


BOOST_AUTO_TEST_SUITE( OptimizerTestSuite )

BOOST_AUTO_TEST_CASE( TestOptimizerConstructor ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4, 5, 2, 1};
    std::vector<double> param2 = {6, 3, 9, 8, 7};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();

    std::vector<double> chebyshevNodes = quadrature.getChebyshevNodes();
    std::vector<double> chebyshevWeights = quadrature.getChebyshevWeights();
    std::vector<std::vector<std::vector<double>>> basisCoefficients = quadrature.getBasisCoefficients();
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis = quadrature.getSplitCompressedBasis();
    std::vector<double> basisIntegrals = quadrature.getBasisIntegrals();

    OptimizerFixture optimizer(chebyshevNodes, chebyshevWeights, basisCoefficients, splitCompressedBasis, basisIntegrals);

    optimizer.testFormJacobian();
} 


BOOST_AUTO_TEST_SUITE_END()