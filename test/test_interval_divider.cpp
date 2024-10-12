#define BOOST_TEST_MODULE IntervalDividerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "interval_divider.hpp"


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
struct IntervalDividerFixture: public IntervalDivider<InputClass>
{
    using InputMethodType = double (InputClass::*)(const double&);
    using InputFunctionType = double(*)(const double&);

    IntervalDividerFixture(int kInput, double lowerBoundInput, double upperBoundInput, InputFunctionType inputFunctionPtr = nullptr, 
    InputMethodType inputMethodPtr = nullptr, InputClass* objectPtr = nullptr) 
    : IntervalDivider<InputClass>(kInput, lowerBoundInput, upperBoundInput, inputFunctionPtr, inputMethodPtr, objectPtr) {};

    ~IntervalDividerFixture() {};

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

    matrix<double> testGetInvertedLegendreMatrix()
    {
        return this->getInvertedLegendreMatrix();
    }

    std::vector<double> testGetLagrangeVector() const
    {
        return this->getLagrangeVector();
    }

    vector<double> testGetAlphaVector() const
    {
        return this->getAlphaVector();
    }

    double testGetMeasure() const
    {
        return this->getMeasure();
    }

    double testGetPrecision() const
    {
        return this->getPrecision();
    }

    inline double testTransformNode(double value) const
    {
        return this->transformNode(value);
    }
};


BOOST_AUTO_TEST_SUITE( DividerTestSuite )

BOOST_AUTO_TEST_CASE( TestDividerValidation ) {
    double lowerBound = 1;
    double upperBound = 2;

    BOOST_CHECK_THROW(IntervalDivider<TestClass> divider(-10, lowerBound, upperBound, testFunction, nullptr, nullptr), std::invalid_argument);
    BOOST_CHECK_NO_THROW(IntervalDivider<TestClass> divider(30, lowerBound, upperBound, testFunction, nullptr, nullptr));
}


BOOST_AUTO_TEST_CASE( TestLegendreMeshEven ) {  
    double lowerBound = -1;
    double upperBound = 1;
    int k = 1;
    std::vector<double> roots = {-0.57735, 0.57735};

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);

    double measure = divider.determineAlphaOnSubinterval();

    std::vector<double> estimatedRoots = divider.testGetLegendreMesh();

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

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);

    double measure = divider.determineAlphaOnSubinterval();

    std::vector<double> estimatedRoots = divider.testGetTransformedMesh();

    BOOST_CHECK_EQUAL(estimatedRoots.size(), 2*k);

    for (size_t i = 0; i < estimatedRoots.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(estimatedRoots[i], roots[i], 1e-6);
    }
}


BOOST_AUTO_TEST_CASE( TestLegendreMatrix ) {
    double lowerBound = 1;
    double upperBound = 2;
    int k = 1;
    matrix<double> legendreMatrix;
    matrix<double> expectedMatrix(2*k, 2*k);

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);

    double measure = divider.determineAlphaOnSubinterval();

    legendreMatrix = divider.testGetLegendreMatrix();

    expectedMatrix(0, 0) = 1;
    expectedMatrix(1, 0) = 1;
    expectedMatrix(0, 1) = divider.testTransformNode(-0.57735);
    expectedMatrix(1, 1) = divider.testTransformNode(0.57735);


    BOOST_CHECK_EQUAL(legendreMatrix.size1(), 2*k);
    BOOST_CHECK_EQUAL(legendreMatrix.size2(), 2*k);

    for (size_t i = 0; i < 2*k; ++i)
    {
        for (size_t j = 0; j < 2*k; ++j)
        {
            BOOST_CHECK_CLOSE_FRACTION(legendreMatrix(i,j), expectedMatrix(i,j), 1e-6);
        }
    }
}


BOOST_AUTO_TEST_CASE( TestLagrangeVectorValues ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 1;

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);
    double measure = divider.determineAlphaOnSubinterval();

    std::vector<double> roots = {-0.57735, 0.57735};
    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    std::vector<double> interpolationPoints = divider.testGetLagrangeVector();

    BOOST_CHECK_EQUAL(interpolationPoints.size(), 2*k);

    for (size_t i = 0; i < interpolationPoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(interpolationPoints[i], expectedInterpolationPoints[i], 1e-6);
    }
}


BOOST_AUTO_TEST_CASE( TestLagrangeVectorValuesWithPassedObject ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 1;
    TestClass testObject;
    TestClass* testObjectPtr = &testObject;

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, nullptr, &TestClass::testMethod, testObjectPtr);
    double measure = divider.determineAlphaOnSubinterval();

    std::vector<double> roots = {-0.57735, 0.57735};
    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    std::vector<double> interpolationPoints = divider.testGetLagrangeVector();

    BOOST_CHECK_EQUAL(interpolationPoints.size(), 2*k);

    for (size_t i = 0; i < interpolationPoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(interpolationPoints[i], expectedInterpolationPoints[i], 1e-6);
    }
}


BOOST_AUTO_TEST_CASE( TestLagrangeVectorValuesNonDefaultDomain ) {
    double lowerBound = 1;
    double upperBound = 2;
    int k = 1;
    std::vector<double> roots = {-0.57735, 0.57735};
    std::transform(roots.begin(), roots.end(), roots.begin(), [&lowerBound, &upperBound](double value)
    { 
        return ((upperBound - lowerBound) * value + (upperBound + lowerBound)) / 2;
    });

    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);
    double measure = divider.determineAlphaOnSubinterval();

    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    std::vector<double> interpolationPoints = divider.testGetLagrangeVector();

    for (size_t i = 0; i < interpolationPoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(interpolationPoints[i], expectedInterpolationPoints[i], 1e-6);
    }
}  


BOOST_AUTO_TEST_CASE( TestMatrixInversion ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 30;
    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);
    double measure = divider.determineAlphaOnSubinterval();

    matrix<double> legendreMatrix = divider.testGetLegendreMatrix();

    BOOST_CHECK_NO_THROW(divider.testGetInvertedLegendreMatrix());

    matrix<double> invertedLegendreMatrix = divider.testGetInvertedLegendreMatrix();

    matrix<double> expectedIdentity = prod(legendreMatrix, invertedLegendreMatrix);

    for (size_t i = 0; i < expectedIdentity.size1(); ++i)
    {
        for (size_t j = 0; j < expectedIdentity.size2(); ++j)
        {   
            if (i == j)
                BOOST_CHECK_CLOSE_FRACTION(expectedIdentity(i, j), 1, 1e-10);
            else
                BOOST_CHECK_SMALL(expectedIdentity(i, j), 1e-10);
        }
    }
}


BOOST_AUTO_TEST_CASE( TestAlphaVectorSolution ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 2;
    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);
    double measure = divider.determineAlphaOnSubinterval();

    vector<double> alphaVector = divider.testGetAlphaVector();

    BOOST_CHECK_EQUAL(alphaVector.size(), 2*k);

    std::vector<double> expectedAlphaVector = {0, 0, 0, 1};

    for (size_t i = 0; i < alphaVector.size(); ++i)
    {
        if (expectedAlphaVector[i] == 0)
            BOOST_CHECK_SMALL(alphaVector[i], 1e-10);
        else
            BOOST_CHECK_CLOSE_FRACTION(alphaVector[i], 1, 1e-10);
    }
}


BOOST_AUTO_TEST_CASE( TestMeasureCalculation ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 30;
    IntervalDividerFixture<TestClass> divider(k, lowerBound, upperBound, testFunction, nullptr, nullptr);
    double measure = divider.determineAlphaOnSubinterval();

    BOOST_CHECK_GT(measure, 0);
}


BOOST_AUTO_TEST_SUITE_END()