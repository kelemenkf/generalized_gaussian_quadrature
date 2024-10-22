#define BOOST_TEST_MODULE FunctionHandlerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "function_handler.hpp"


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


BOOST_AUTO_TEST_CASE( TestFunctionProcessing0Param ) {;
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
    testValue = handler.callFunction(5);

    BOOST_CHECK_EQUAL(testValue, 56.0);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessingInvalidParam ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(testFunction3Params, param1, param2);

    BOOST_CHECK_THROW(handler.callFunction(5), std::runtime_error);
}


BOOST_AUTO_TEST_CASE( TestParameterCombinationBuilder ) {
    std::vector<double> param1 = {2, 3};
    std::vector<double> param2 = {1, 25};
    std::vector<double> param3 = {1, 4};

    std::vector<std::vector<double>> expectedCombos1Param = {{2},{3}};

    std::vector<std::vector<double>> expectedCombos2Param = {
        {2, 1},
        {2, 25},
        {3, 1},
        {3, 25}
    };

    std::vector<std::vector<double>> expectedCombos3Param = {
        {2, 1, 1},
        {2, 1, 4},
        {2, 25, 1},
        {2, 25, 4},
        {3, 1, 1},
        {3, 1, 4},
        {3, 25, 1},
        {3, 25, 4}
    };

    FunctionHandler<> handler0(testFunction);
    FunctionHandler<std::vector<double>> handler1(testFunction1Param, param1);
    FunctionHandler<std::vector<double>, std::vector<double>> handler2(testFunction2Param, param1, param2);
    FunctionHandler<std::vector<double>, std::vector<double>, std::vector<double>> handler3(testFunction2Param, param1, param2, param3);
   
    std::vector<std::vector<double>> combo1 = handler1.getParameterCombinations();
    std::vector<std::vector<double>> combo2 = handler2.getParameterCombinations();
    std::vector<std::vector<double>> combo3 = handler3.getParameterCombinations();
    BOOST_CHECK_EQUAL(combo1.size(), 2);
    BOOST_CHECK_EQUAL(combo2.size(), 4);
    BOOST_CHECK_EQUAL(combo3.size(), 8);


    for (size_t i = 0; i < combo1.size(); ++i)
    {
        BOOST_CHECK_EQUAL_COLLECTIONS(combo1[i].begin(), combo1[i].end(), expectedCombos1Param[i].begin(), expectedCombos1Param[i].end());
    }
    for (size_t i = 0; i < combo2.size(); ++i)
    {
        BOOST_CHECK_EQUAL_COLLECTIONS(combo2[i].begin(), combo2[i].end(), expectedCombos2Param[i].begin(), expectedCombos2Param[i].end());
    }
    for (size_t i = 0; i < combo2.size(); ++i)
    {
        BOOST_CHECK_EQUAL_COLLECTIONS(combo3[i].begin(), combo3[i].end(), expectedCombos3Param[i].begin(), expectedCombos3Param[i].end());
    }
}   


BOOST_AUTO_TEST_CASE( TestCombinationIndexIncrementer ) {
    FunctionHandler<> handler(testFunction);

    handler.incrementCombinationIndex();

    size_t index = handler.getCombinationIndex();
    size_t expectedIndex = 1;

    BOOST_CHECK_EQUAL(index, expectedIndex);
}


BOOST_AUTO_TEST_CASE( TestResetCombinationIndex ) {
    FunctionHandler<> handler(testFunction);

    handler.incrementCombinationIndex();
    handler.resetCombinationIndex();

    size_t index = handler.getCombinationIndex();
    size_t expectedIndex = 0;

    BOOST_CHECK_EQUAL(index, expectedIndex);
}


BOOST_AUTO_TEST_SUITE_END()