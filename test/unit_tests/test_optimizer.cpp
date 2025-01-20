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

    MatrixXd testShermanMorrisonWoodburry(const MatrixXd& input, int j)
    {
       VectorXd u_k = input.col(j); 

       MatrixXd rank1 =  u_k * u_k.transpose();

       return input - rank1; 
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
    
    std::vector<std::vector<double>> Jacobian = optimizer.getJacobian();

    BOOST_CHECK_EQUAL(Jacobian.size(), splitCompressedBasis.size());
    BOOST_CHECK_EQUAL(Jacobian[0].size(), 2 * chebyshevNodes.size());
}


// BOOST_AUTO_TEST_CASE( TestShermannWoodburryMorrision ) {
//     MatrixXd A {
//         {1, 3},
//         {2, 4}
//     };

//     MatrixXd expected {
//         {0, 1},
//         {0, 0}
//     };

//     double lowerBound = 0;
//     double upperBound = 2; 
//     int k = 30;
//     std::vector<double> param1 = {5, 4};
//     std::vector<double> param2 = {6, 3};
//     FunctionHandler<std::vector<double>, std::vector<double>> handlerTest3Param(testFunction2ParamPC, param1, param2);
    
//     QuadratureRule quadrature(lowerBound, upperBound, handlerTest3Param);
//     quadrature.calculateQuadratureNodes();
//     quadrature.compressFunctionSpace();
//     quadrature.obtainBasisCoefficients();
//     quadrature.optimizeQuadrature();

//     std::vector<double> chebyshevNodes = quadrature.getChebyshevNodes();
//     std::vector<double> chebyshevWeights = quadrature.getChebyshevWeights();
//     std::vector<std::vector<std::vector<double>>> basisCoefficients = quadrature.getBasisCoefficients();
//     std::vector<std::vector<std::vector<double>>> splitCompressedBasis = quadrature.getSplitCompressedBasis();
//     std::vector<double> basisIntegrals = quadrature.getBasisIntegrals();
//     std::vector<double> endpoints = quadrature.getConsolidatedEndpoints();
//     std::vector<std::vector<double>> splitNodes = quadrature.getSplitNodes();

//     OptimizerFixture optimizer(chebyshevNodes, chebyshevWeights, basisCoefficients, splitCompressedBasis, basisIntegrals, endpoints, splitNodes);

//     MatrixXd calculated = optimizer.testShermanMorrisonWoodburry(A, 0);     

//     for (size_t row = 0; row < expected.rows(); ++row)
//     {
//         for (size_t col = 0; col < expected.cols(); ++col)
//         {
//             BOOST_CHECK_EQUAL(calculated(row, col), expected(row, col));
//         }
//     }

//     std::cout << calculated << std::endl;
// }


BOOST_AUTO_TEST_CASE( TestCalculateStepDirections ) {
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

    std::vector<VectorXd> stepDirections = optimizer.getStepDirections();

    for (size_t i = 0; i < stepDirections.size(); ++i)
    {
        std::cout << stepDirections[i] << std::endl;
    }
}


BOOST_AUTO_TEST_SUITE_END()