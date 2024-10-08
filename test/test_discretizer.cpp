#define BOOST_TEST_MODULE DiscretizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "discretizer.hpp"


class TestClass
{
public:
    double testMethod(const std::vector<double>& input)
    {
        std::cout << "Test method ran" << std::endl;
        return input[0];
    }
};


double testFunction(const std::vector<double>& input)
{
    std::cout << "Test function ran" << std::endl;
    return input[0];
}


BOOST_AUTO_TEST_SUITE( DiscretizerTestSuite )

BOOST_AUTO_TEST_CASE( TestDiscretizerValidation ) {
    double lowerBound = 1;
    double upperBound = 2;
    
    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(-10, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(30, -0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_NO_THROW(Discretizer<TestClass> discretizer(30, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr));
}

BOOST_AUTO_TEST_SUITE_END()