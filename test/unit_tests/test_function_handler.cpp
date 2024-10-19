#define BOOST_TEST_MODULE FunctionHandlerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "function_handler.hpp"


using ParameterSpace = std::variant
<
    std::tuple<double>,
    std::tuple<double, double>,
    std::tuple<double, double, double>
>;


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


double testFunction1Param(const double& x, const double& param1)
{
    return param1 * x;
}


double testFunction2Param(const double& x, const double& param1, const double& param2)
{
    return param1 * x + param2;
}


double testFunction3Params(const double& x, const double& param1, const double& param2, const double& param3)
{
    return param1 * (x*x) + param2 * x + param3;
}


BOOST_AUTO_TEST_SUITE( FunctionHandlerTestSuite )


BOOST_AUTO_TEST_CASE( TestParamSpaceConstruction ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};

    FunctionHandler<std::vector<double>> handler(testFunction, param1);

    std::vector<std::vector<double>> paramSpace = handler.getParamSpace();

    BOOST_CHECK_EQUAL(paramSpace.size(), 1);
    BOOST_CHECK_EQUAL(paramSpace[0].size(), 6);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessing0Param ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};

    FunctionHandler<> handler(testFunction);

    double testValue;
    testValue = handler.callFunction(5);

    BOOST_CHECK_EQUAL(testValue, 5.0);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessing3Param ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};
    std::vector<double> param3 = {1, 4, 6, 43, 56, 90};

    FunctionHandler<std::vector<double>, std::vector<double>, std::vector<double>> handler(testFunction3Params, param1, param2, param3);

    double testValue;
    testValue = handler.callFunction(5, 2, 4, -6);

    BOOST_CHECK_EQUAL(testValue, 64.0);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessingInvalidParam ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(testFunction3Params, param1, param2);

    BOOST_CHECK_THROW(handler.callFunction(5, 2, 4, -6), std::runtime_error);
}


BOOST_AUTO_TEST_CASE( TestParameterCombinationBuilder ) {
    std::vector<double> param1 = {2, 3};
    std::vector<double> param2 = {1, 25};
    std::vector<double> param3 = {1, 4};

    std::vector<ParameterSpace> expectedCombos1Param = {
        ParameterSpace(std::make_tuple(2)),
        ParameterSpace(std::make_tuple(3))
    };

    std::vector<ParameterSpace> expectedCombos2Param = {
        ParameterSpace(std::make_tuple(2, 1)),
        ParameterSpace(std::make_tuple(2, 25)),
        ParameterSpace(std::make_tuple(3, 1)),
        ParameterSpace(std::make_tuple(3, 25))
    };

    std::vector<ParameterSpace> expectedCombos3Param = {
        ParameterSpace(std::make_tuple(2, 1, 1)),
        ParameterSpace(std::make_tuple(2, 1, 4)),
        ParameterSpace(std::make_tuple(2, 25, 1)),
        ParameterSpace(std::make_tuple(2, 25, 4)),
        ParameterSpace(std::make_tuple(3, 1, 1)),
        ParameterSpace(std::make_tuple(3, 1, 4)),
        ParameterSpace(std::make_tuple(3, 25, 1)),
        ParameterSpace(std::make_tuple(3, 25, 4))
    };

    FunctionHandler<> handler0(testFunction);
    FunctionHandler<std::vector<double>> handler1(testFunction1Param, param1);
    FunctionHandler<std::vector<double>, std::vector<double>> handler2(testFunction2Param, param1, param2);
    FunctionHandler<std::vector<double>, std::vector<double>, std::vector<double>> handler3(testFunction2Param, param1, param2, param3);
   

    handler1.buildParameterCombinations();
    handler2.buildParameterCombinations();
    handler3.buildParameterCombinations();
    std::vector<ParameterSpace> combo1 = handler1.getParameterCombinations();
    std::vector<ParameterSpace> combo2 = handler2.getParameterCombinations();
    std::vector<ParameterSpace> combo3 = handler3.getParameterCombinations();
    BOOST_CHECK_EQUAL(combo1.size(), 2);
    BOOST_CHECK_EQUAL(combo2.size(), 4);
    BOOST_CHECK_EQUAL(combo3.size(), 8);
}   


BOOST_AUTO_TEST_SUITE_END()