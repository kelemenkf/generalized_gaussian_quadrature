#define BOOST_TEST_MODULE DiscretizerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "discretizer.hpp"


class TestClass
{
public:
    double testMethod(const double& x)
    {
        return (5*std::pow(x,3) - 3*x) / 2;
    }
};


double testFunction(const double& x)
{
    return (5*std::pow(x,3) - 3*x) / 2;
}


template<typename InputClass>
struct DiscretizerFixture: public Discretizer<InputClass>
{
    using InputMethodType = double (InputClass::*)(const double&);
    using InputFunctionType = double(*)(const double&);

    DiscretizerFixture(int kInput, double precisionInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : Discretizer<InputClass>(kInput, precisionInput, lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr) {};

    ~DiscretizerFixture() {};

    void testDiscretizationRoutione() 
    {
        this->discretizationRoutine();
    }

    double testGetPrecision() const
    {
        return this->getPrecision();
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

BOOST_AUTO_TEST_SUITE_END()