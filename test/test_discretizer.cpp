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

    std::vector<double> testGetTransformedMesh() const
    {
        return this->getTransformedMesh();
    }

    matrix<double> testGetLegendreMatrix()
    {
        return this->getLegendreMatrix();
    }

    inline double transformNode(double value) const
    {
        return this->transformNode(value);
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
    double lowerBound = -1;
    double upperBound = 1;
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


BOOST_AUTO_TEST_CASE( TestLegendreMeshTransformed ) {  
    double lowerBound = 1;
    double upperBound = 2;
    int k = 1;
    std::vector<double> roots = {-0.57735, 0.57735};
    std::transform(roots.begin(), roots.end(), roots.begin(), [&lowerBound, &upperBound](double value)
    { 
        return ((upperBound - lowerBound) * value + (upperBound + lowerBound)) / 2;
    });

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    std::vector<double> estimatedRoots = discretizer.testGetTransformedMesh();

    BOOST_CHECK_EQUAL(estimatedRoots.size(), 2*k);

    for (size_t i = 0; i < estimatedRoots.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(estimatedRoots[i], roots[i], 0.0001);
    }
}


BOOST_AUTO_TEST_CASE( TestLegendreMatrix ) {
    double lowerBound = 1;
    double upperBound = 2;
    int k = 1;
    matrix<double> legendreMatrix;
    matrix<double> expectedMatrix(2*k, 2*k);

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    legendreMatrix = discretizer.testGetLegendreMatrix();

    expectedMatrix(0, 0) = 1;
    expectedMatrix(0, 1) = 1;
    expectedMatrix(1, 0) = -0.57735;
    expectedMatrix(1, 1) = 0.57735;


    BOOST_CHECK_EQUAL(legendreMatrix.size1(), 2*k);
    BOOST_CHECK_EQUAL(legendreMatrix.size2(), 2*k);

    for (size_t i = 0; i < 2*k; ++i)
    {
        for (size_t j = 0; j < 2*k; ++j)
        {
            BOOST_CHECK_EQUAL(legendreMatrix(i,j), expectedMatrix(i,j));
        }
    }
}


BOOST_AUTO_TEST_SUITE_END()