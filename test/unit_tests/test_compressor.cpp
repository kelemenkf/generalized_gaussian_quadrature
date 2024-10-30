#define BOOST_TEST_MODULE CompressorTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "compressor.hpp"
#include "ggq.hpp"


double testFunction2ParamPC(const double& x, const double& param1, const double& param2)
{
    if (x <= 1)
    {
        return param1 * x * x;
    }
    else 
    {
        return param2 - x;
    }
}


bool isOrthonormal(const Eigen::MatrixXd& matrix, double tolerance = 1e-9) {
    if (matrix.rows() != matrix.cols()) {
        return false;
    }

    Eigen::MatrixXd identity_check = matrix.transpose() * matrix;

    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols());
    return (identity_check.isApprox(identity, tolerance));
}


bool isOrthonormalWithWeight(const Eigen::MatrixXd& matrix, const std::vector<double>& weights, double tolerance = 1e-9)
{
    if (matrix.rows() != matrix.cols()) {
        return false;
    }

    const Eigen::Map<const Eigen::VectorXd> weightsEigen(weights.data(), weights.size());

    Eigen::MatrixXd identity_check = matrix.transpose() * weightsEigen.asDiagonal() * matrix;

    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols());
    return (identity_check.isApprox(identity, tolerance));
}


double lowerBound = 0;
double upperBound = 2; 
std::vector<double> param1 = {5, 4};
std::vector<double> param2 = {6, 3};
FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);

QuadratureRule<FunctionHandler<std::vector<double>, std::vector<double>>> quadratureObject(lowerBound, upperBound, handlerPiecewiseSmooth);
QuadratureRule<FunctionHandler<std::vector<double>, std::vector<double>>>* quadrature = &quadratureObject;

template<typename T>
struct CompressorFixture: public Compressor<T>
{
    CompressorFixture(const T* quadratureInput, double quadraturePrecisionInput = 1e-3) : 
    Compressor<T>(quadratureInput, quadraturePrecisionInput) {}

    ~CompressorFixture() {}

    void testConstructA()
    {
        this->constructA();
    }
}; 


BOOST_AUTO_TEST_SUITE( CompressorTestSuite )


BOOST_AUTO_TEST_CASE( TestCompressorConstructor ) {
    CompressorFixture compressor(quadrature);
}


BOOST_AUTO_TEST_CASE( TestConstructA ) {
    quadrature->calculateQuadratureNodes();
    CompressorFixture compressor(quadrature);

    MatrixXd A = compressor.getA();

    BOOST_CHECK_EQUAL(A.rows(), 60);
    BOOST_CHECK_EQUAL(A.cols(), 4);
}


BOOST_AUTO_TEST_CASE( TestQRDecompositon ) {
    CompressorFixture compressor(quadrature);

    MatrixXd U = compressor.getU();

    BOOST_CHECK_EQUAL(isOrthonormal(U), true);
}


BOOST_AUTO_TEST_CASE( TestUScaling ) {
    CompressorFixture compressor(quadrature);

    MatrixXd scaledU = compressor.getScaledU();
    MatrixXd originalU = compressor.getU();
    std::vector<double> weights = quadrature->getWeights();

    MatrixXd U = scaledU;
    for (size_t row = 0; row < U.rows(); ++row)
    {
        U.row(row) *= sqrt(weights[row]);
    }

    BOOST_CHECK_EQUAL(U.rows(), originalU.rows());
    BOOST_CHECK_EQUAL(U.cols(), originalU.cols());


    for (size_t row = 0; row < U.rows(); ++row)
    {
        for (size_t column = 0; column < U.cols(); ++column)
        {
            BOOST_CHECK_CLOSE_FRACTION(originalU(row, column), U(row, column), 1e-6);
        }
    }

    BOOST_CHECK_EQUAL(isOrthonormalWithWeight(scaledU, weights), true);
}


BOOST_AUTO_TEST_CASE( TestLegendreWeights ) {
    CompressorFixture compressor(quadrature);

    std::vector<double> weights = quadrature->getWeights();

    double weightSumFirstHalf = std::accumulate(weights.begin(), weights.begin() + (weights.size() / 2), 0.0);

    double weightSumSecondHalf = std::accumulate(weights.begin() + (weights.size() / 2), weights.end(), 0.0);


    BOOST_CHECK_CLOSE_FRACTION(weightSumFirstHalf, 1.0, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(weightSumSecondHalf, 1.0, 1e-6); 
}


BOOST_AUTO_TEST_CASE( TestGetNormalizingFactors ) {
    CompressorFixture compressor(quadrature);

    std::vector<double> normalizingFactors = compressor.getNormalizingFactors();
    MatrixXd A = compressor.getA();

    BOOST_CHECK_EQUAL(normalizingFactors.size(), A.cols());

    for (size_t i = 0; i < normalizingFactors.size(); ++i)
    {
        BOOST_CHECK_GT(normalizingFactors[i], 0);
    }
}


BOOST_AUTO_TEST_CASE( TestOrthonomalBasis ) {
    CompressorFixture compressor(quadrature);

    std::vector<double> normalizingFactors = compressor.getNormalizingFactors();
    std::vector<std::vector<double>> compressedBasis = compressor.getCompressedBasis();

    BOOST_CHECK_EQUAL(compressedBasis.size(), 3);
}



BOOST_AUTO_TEST_CASE( TestRVector ) {
    CompressorFixture compressor(quadrature);

    std::vector<double> rVector = compressor.getRVector();

    displayVector(rVector);

    BOOST_CHECK_EQUAL(rVector.size(), 3);
}


BOOST_AUTO_TEST_CASE( TestConstructB ) {
    CompressorFixture compressor(quadrature);

    MatrixXd B = compressor.getB();

    BOOST_CHECK_EQUAL(B.rows(), 3);
    BOOST_CHECK_EQUAL(B.cols(), 60);
}


BOOST_AUTO_TEST_SUITE_END()