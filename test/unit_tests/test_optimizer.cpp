#define BOOST_TEST_MODULE OptimizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/algorithm/cxx11/is_sorted.hpp>
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

    MatrixXd testShermanMorrisonWoodburry(const MatrixXd& input, int j)
    {
       VectorXd u_k = input.col(j); 

       MatrixXd rank1 =  u_k * u_k.transpose();

       return input - rank1; 
    }
};

std::vector<double> param1 = {5, 4, 5, 2, 1};
std::vector<double> param2 = {6, 3, 9, 8, 7}; 
FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, true, param1, param2);


BOOST_AUTO_TEST_SUITE( OptimizerTestSuite )

BOOST_AUTO_TEST_CASE( TestOptimizerConstructor ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4, 5, 2, 1};
    std::vector<double> param2 = {6, 3, 9, 8, 7};
    
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
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
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

    std::map<int, std::vector<std::pair<double, double>>> intervalNodeMap = optimizer.getIntervalChebyshevNodesMap();

    int numberOfNodes = 0;

    for (auto element = intervalNodeMap.begin(); element != intervalNodeMap.end(); ++element)
    {
        std::cout << element->first << std::endl; 
        numberOfNodes += (element->second).size();
    }

    BOOST_CHECK_EQUAL(numberOfNodes, chebyshevNodes.size());
}


BOOST_AUTO_TEST_CASE( TestFormJacobian ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
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
    
    std::vector<std::vector<double>> Jacobian = optimizer.getJacobian();

    BOOST_CHECK_EQUAL(Jacobian.size(), splitCompressedBasis.size());
    BOOST_CHECK_EQUAL(Jacobian[0].size(), 2 * chebyshevNodes.size());
}


BOOST_AUTO_TEST_CASE( TestCalculateStepDirections ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
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

    std::vector<VectorXd> stepDirections = optimizer.getStepDirections();

    int n = chebyshevNodes.size(); 

    BOOST_CHECK_EQUAL(stepDirections.size(), n);

    std::vector<double> stepDirectionNorms = optimizer.getStepDirectionNorms();

    BOOST_CHECK_EQUAL(stepDirectionNorms.size(), n);
}


BOOST_AUTO_TEST_CASE( TestStepDirectionNormReordering ) {
    double lowerBound = 0;
    double upperBound = 2; 
    int k = 30;
    std::vector<double> param1 = {5, 4};
    std::vector<double> param2 = {6, 3};
    
    QuadratureRule quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);
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

    std::vector<std::pair<double, int>> sortedStepDirectionNorms = optimizer.getSortedStepDirectionNorms();

    std::vector<double> sortedNorms; 

    for (auto element = sortedStepDirectionNorms.begin(); element != sortedStepDirectionNorms.end(); ++element)
    {
        sortedNorms.push_back(element->first);
    }

    BOOST_CHECK_EQUAL(boost::algorithm::is_sorted(sortedNorms), true);
}


BOOST_AUTO_TEST_SUITE_END()