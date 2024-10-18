#define BOOST_TEST_MODULE IntervalDividerTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "interval_divider.hpp"
#include "function_handler.hpp"


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


FunctionHandler<> handler(testFunction);


template<typename T>
struct IntervalDividerFixture: public IntervalDivider<T>
{
    IntervalDividerFixture(int kInput, double lowerBoundInput, double upperBoundInput, const T& handlerInput) 
    : IntervalDivider<T>(kInput, lowerBoundInput, upperBoundInput, handlerInput) {};

    ~IntervalDividerFixture() {};

    std::vector<double> testGetLegendreMesh() const
    {
        return this->getLegendreMesh();
    }

    std::vector<double> testGetTransformedMesh() const
    {
        return this->getTransformedMesh();
    }

    matrix<double> testGetLegendreMatrix() const
    {
        return this->getLegendreMatrix();
    }

    matrix<double> testGetInvertedLegendreMatrix() const
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

    inline double testTransformNode(double value) const
    {
        return this->transformNode(value);
    }
};


BOOST_AUTO_TEST_SUITE( DividerTestSuite )

BOOST_AUTO_TEST_CASE( TestDividerValidation ) {
    double lowerBound = 1;
    double upperBound = 2;

    BOOST_CHECK_NO_THROW(IntervalDivider divider(30, lowerBound, upperBound, handler));
}


BOOST_AUTO_TEST_CASE( TestLegendreMeshEven ) {  
    double lowerBound = -1;
    double upperBound = 1;
    int k = 1;
    std::vector<double> roots = {-0.57735, 0.57735};

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);

    divider.processInterval();

    std::vector<double> estimatedRoots = divider.testGetLegendreMesh();

    BOOST_CHECK_EQUAL(estimatedRoots.size(), 2*k);

    for (size_t i = 0; i < estimatedRoots.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(estimatedRoots[i], roots[i], 1e-6);
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

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);

    divider.processInterval();

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

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);

    divider.processInterval();

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

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

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

    auto method = std::bind(&TestClass::testMethod, testObjectPtr, std::placeholders::_1);
    std::function<double(const double&)> methodPtr = method;
    FunctionHandler<> handler(methodPtr);

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

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

    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

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
    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

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
    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

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
    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

    double measure = divider.getMeasure();

    BOOST_CHECK_GT(measure, 0);
}


BOOST_AUTO_TEST_CASE( TestMeasureCalculationNonStandardInterval ) {
    double lowerBound = 1;
    double upperBound = 2;
    int k = 30;
    IntervalDividerFixture divider(k, lowerBound, upperBound, handler);
    divider.processInterval();

    double measure = divider.getMeasure();

    BOOST_CHECK_GT(measure, 0);
}


BOOST_AUTO_TEST_SUITE_END()