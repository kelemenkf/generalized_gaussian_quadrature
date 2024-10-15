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
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(testFunction, param1, param2);

    std::vector<std::vector<double>> paramSpace = handler.getParamSpace();

    BOOST_CHECK_EQUAL(paramSpace.size(), 2);
    BOOST_CHECK_EQUAL(paramSpace[0].size(), 6);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessing1Param ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(testFunction, param1, param2);

    double testValue;
    testValue = handler.callFunction(5);

    BOOST_CHECK_EQUAL(testValue, 5.0);
}


BOOST_AUTO_TEST_CASE( TestFunctionProcessing3Param ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(testFunction3Params, param1, param2);

    double testValue;
    testValue = handler.callFunction(5, 2, 4, -6);

    BOOST_CHECK_EQUAL(testValue, 64.0);
}


BOOST_AUTO_TEST_SUITE_END()