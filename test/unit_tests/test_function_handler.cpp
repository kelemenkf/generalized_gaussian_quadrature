#define BOOST_TEST_MODULE FunctionHandlerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "function_handler.hpp"


BOOST_AUTO_TEST_SUITE( FunctionHandlerTestSuite )

BOOST_AUTO_TEST_CASE( TestParamSpaceConstruction ) {
    std::vector<double> param1 = {2, 3, 5, 6, 7, 0};
    std::vector<double> param2 = {1, 25, 7, 3, 6, 2};

    FunctionHandler<std::vector<double>, std::vector<double>> handler(param1, param2);

    std::vector<std::vector<double>> paramSpace = handler.getParamSpace();

    BOOST_CHECK_EQUAL(paramSpace.size(), 2);
    BOOST_CHECK_EQUAL(paramSpace[0].size(), 6);
}

BOOST_AUTO_TEST_SUITE_END()