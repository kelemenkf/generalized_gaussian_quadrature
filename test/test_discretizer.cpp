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
        BOOST_CHECK_CLOSE_FRACTION(estimatedRoots[i], roots[i], 1e-6);
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
    expectedMatrix(1, 0) = 1;
    expectedMatrix(0, 1) = discretizer.testTransformNode(-0.57735);
    expectedMatrix(1, 1) = discretizer.testTransformNode(0.57735);


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

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    std::vector<double> roots = {-0.57735, 0.57735};
    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    discretizer.testDiscretizationRoutione();
    std::vector<double> interpolationPoints = discretizer.testGetLagrangeVector();

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

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, nullptr, &TestClass::testMethod, testObjectPtr);
    std::vector<double> roots = {-0.57735, 0.57735};
    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    discretizer.testDiscretizationRoutione();
    std::vector<double> interpolationPoints = discretizer.testGetLagrangeVector();

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

    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    std::vector<double> expectedInterpolationPoints = {testFunction(roots[0]), testFunction(roots[1])};
    discretizer.testDiscretizationRoutione();
    std::vector<double> interpolationPoints = discretizer.testGetLagrangeVector();

    for (size_t i = 0; i < interpolationPoints.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(interpolationPoints[i], expectedInterpolationPoints[i], 1e-6);
    }
}  


BOOST_AUTO_TEST_CASE( TestMatrixInversion ) {
    double lowerBound = -1;
    double upperBound = 1;
    int k = 30;
    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);
    matrix<double> legendreMatrix = discretizer.testGetLegendreMatrix();

    BOOST_CHECK_NO_THROW(discretizer.testGetInvertedLegendreMatrix());

    matrix<double> invertedLegendreMatrix = discretizer.testGetInvertedLegendreMatrix();

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
    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);

    discretizer.discretizationRoutine();

    vector<double> alphaVector = discretizer.testGetAlphaVector();

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
    DiscretizerFixture<TestClass> discretizer(k, 0.01, lowerBound, upperBound, testFunction, nullptr, nullptr);

    discretizer.discretizationRoutine();

    double measure = discretizer.testGetMeasure();

    std::cout << discretizer.testGetPrecision() << std::endl;

    BOOST_CHECK_GT(measure, 0);
}


BOOST_AUTO_TEST_SUITE_END()