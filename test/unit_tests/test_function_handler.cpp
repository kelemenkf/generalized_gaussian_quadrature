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


BOOST_AUTO_TEST_CASE( TestFunctionProcessing1Param ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};

    FunctionHandler<std::vector<double>> handler(testFunction, param1);

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

    double testValue;
    BOOST_CHECK_THROW(testValue = handler.callFunction(5, 2, 4, -6), std::runtime_error);
}


BOOST_AUTO_TEST_SUITE_END()