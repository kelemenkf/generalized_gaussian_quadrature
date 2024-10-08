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


template<typename InputClass>
struct DiscretizerFixture: public Discretizer<InputClass>
{
    using InputMethodType = double (InputClass::*)(const std::vector<double>&);
    using InputFunctionType = double(*)(const std::vector<double>&);

    DiscretizerFixture(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : Discretizer<InputClass>(kInput, precisionInput, lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr) {};

    ~DiscretizerFixture() {};

    std::vector<double> testGetLegendreMesh() const
    {
        return this->getLegendreMesh();
    }
};


BOOST_AUTO_TEST_SUITE( DiscretizerTestSuite )

BOOST_AUTO_TEST_CASE( TestDiscretizerValidation ) {
    double lowerBound = 1;
    double upperBound = 2;

    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(-10, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_THROW(Discretizer<TestClass> discretizer(30, -0.01, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_NO_THROW(Discretizer<TestClass> discretizer(30, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr));
}


BOOST_AUTO_TEST_CASE( TestLegendreMeshEven ) {  
    double lowerBound = 1;
    double upperBound = 2;
    int k = 1;
    std::vector<double> roots = {-0.57735, 0.57735};

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    std::vector<double> estimatedRoots = discretizer.testGetLegendreMesh();

    BOOST_CHECK_EQUAL(estimatedRoots.size(), 2*k);

    for (size_t i = 0; i < estimatedRoots.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(estimatedRoots[i], roots[i], 0.0001);
    }
}

BOOST_AUTO_TEST_SUITE_END()