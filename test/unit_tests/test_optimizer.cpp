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


double testFunction3Param(const double& x, const double& param1, const double& param2, const double& param3)
{
   if (x < 10)
   {
        return param1 * x * x;
   } 
   else if (x >= 10 && x < 20)
   {
        return param2 - x; 
   }
   else 
   {
        return - sin(param3 * x);
   }
}


struct OptimizerFixture : public Optimizer
{
    OptimizerFixture(const std::vector<double>& inputChebyshevNodes, const std::vector<double>& inputChebyshevWeights, 
    const std::vector<std::vector<std::vector<double>>>& inputBasisCoefficients, const std::vector<std::vector<std::vector<double>>>& 
    inputSplitCompressedBasis, const std::vector<double>& inputBasisIntegrals, const std::vector<double>& inputEndpoints,
    const std::vector<std::vector<double>> splitNodes)
    : Optimizer(inputChebyshevNodes, inputChebyshevWeights, inputBasisCoefficients, inputSplitCompressedBasis, inputBasisIntegrals,
    inputEndpoints, splitNodes) 
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
    std::vector<double> endpoints = quadrature.getConsolidatedEndpoints();
    std::vector<std::vector<double>> splitNodes = quadrature.getSplitNodes();

    displayVector(chebyshevNodes);
} 


BOOST_AUTO_TEST_CASE( TestOptimizerIntervalChyebyshevNodeMap ) {    
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerTest3Param(testFunction2ParamPC, param1, param2);
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerTest3Param);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();
    quadrature.optimizeQuadrature();

    std::vector<double> chebyshevNodes = quadrature.getChebyshevNodes();
    std::vector<double> chebyshevWeights = quadrature.getChebyshevWeights();
    std::vector<std::vector<std::vector<double>>> basisCoefficients = quadrature.getBasisCoefficients();
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis = quadrature.getSplitCompressedBasis();
    std::vector<double> basisIntegrals = quadrature.getBasisIntegrals();
    std::vector<double> endpoints = quadrature.getConsolidatedEndpoints();
    std::vector<std::vector<double>> splitNodes = quadrature.getSplitNodes();

    OptimizerFixture optimizer(chebyshevNodes, chebyshevWeights, basisCoefficients, splitCompressedBasis, basisIntegrals, endpoints, splitNodes);

    //std::map<int, std::vector<double>> intervalNodeMap = optimizer.getIntervalChebyshevNodesMap();

    //for (auto element = intervalNodeMap.begin(); element != intervalNodeMap.end(); ++element)
    //{
    //    std::cout << element->first << std::endl; 
    //   displayVector(element->second);
    //}
}


BOOST_AUTO_TEST_CASE( TestFormJacobian ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    FunctionHandler<std::vector<double>, std::vector<double>> handlerTest3Param(testFunction2ParamPC, param1, param2);
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerTest3Param);
    quadrature.calculateQuadratureNodes();
    quadrature.compressFunctionSpace();
    quadrature.obtainBasisCoefficients();
    quadrature.optimizeQuadrature();

    std::vector<double> chebyshevNodes = quadrature.getChebyshevNodes();
    std::vector<double> chebyshevWeights = quadrature.getChebyshevWeights();
    std::vector<std::vector<std::vector<double>>> basisCoefficients = quadrature.getBasisCoefficients();
    std::vector<std::vector<std::vector<double>>> splitCompressedBasis = quadrature.getSplitCompressedBasis();
    std::vector<double> basisIntegrals = quadrature.getBasisIntegrals();
    std::vector<double> endpoints = quadrature.getConsolidatedEndpoints();
    std::vector<std::vector<double>> splitNodes = quadrature.getSplitNodes();

    OptimizerFixture optimizer(chebyshevNodes, chebyshevWeights, basisCoefficients, splitCompressedBasis, basisIntegrals, endpoints, splitNodes);
}


BOOST_AUTO_TEST_SUITE_END()