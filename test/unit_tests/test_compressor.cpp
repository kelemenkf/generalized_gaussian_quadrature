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


MatrixXd reorderMatrix(const MatrixXd& B, const std::vector<size_t>& indices) {
    MatrixXd reorderedB(B.rows(), B.cols());

    for (size_t i = 0; i < indices.size(); ++i) {
        size_t originalIndex = indices[i]; 
        reorderedB.col(originalIndex) = B.col(i); 
    }

    return reorderedB;
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

    std::tuple<MatrixXd, MatrixXd, std::vector<size_t>> testPivotedGramSchmidt(MatrixXd& inputMatrix)
    {   
        return this->doublePivotedGramSchmidt(inputMatrix);
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

    VectorXd rVector = compressor.getRVector();

    BOOST_CHECK_EQUAL(rVector.size(), 3);
}


BOOST_AUTO_TEST_CASE( TestConstructB ) {
    CompressorFixture compressor(quadrature);

    MatrixXd B = compressor.getB();

    BOOST_CHECK_EQUAL(B.rows(), 3);
    BOOST_CHECK_EQUAL(B.cols(), 60);
}


BOOST_AUTO_TEST_CASE( TestPivotedGramSchmidt ) {
    CompressorFixture compressor(quadrature);

    MatrixXd B{
        {1, 1, 0, 3},
        {1, 2, 0, 4},
        {0, 0, 1, 2}
    };
    std::tuple<MatrixXd, MatrixXd, std::vector<size_t>> result = compressor.testPivotedGramSchmidt(B);

    MatrixXd Q = std::get<0>(result);
    MatrixXd R = std::get<1>(result);
    std::vector<size_t> perm = std::get<2>(result);

    MatrixXd computedB = Q * R;

    MatrixXd reorderedB = reorderMatrix(computedB, perm);

    if (B.cols() == reorderedB.cols() && B.rows() == reorderedB.rows())
    {
        for (size_t i = 0; i < B.rows(); ++i)
        {
            for (size_t j = 0; j < B.cols(); ++j)
            {
                if (B(i,j) == 0)
                {
                    BOOST_CHECK_SMALL(reorderedB(i,j), 1e-6);
                }
                else
                {
                    BOOST_CHECK_CLOSE_FRACTION(B(i,j), reorderedB(i,j), 1e-6);
                }
            }      
        }
    }
    else
    {
        std::cout << "Sizes don't match" << std::endl;
    }
}


 
BOOST_AUTO_TEST_CASE( TestPivotedGramSchmidtWithB ) {
    CompressorFixture compressor(quadrature);

    MatrixXd B = compressor.getB();
    std::tuple<MatrixXd, MatrixXd, std::vector<size_t>> result = compressor.testPivotedGramSchmidt(B);

    MatrixXd Q = std::get<0>(result);
    MatrixXd R = std::get<1>(result);
    std::vector<size_t> perm = std::get<2>(result);

    BOOST_CHECK_EQUAL(isOrthonormal(Q), true);

    MatrixXd computedB = Q * R;

    MatrixXd reorderedB = reorderMatrix(computedB, perm);

    if (B.cols() == reorderedB.cols() && B.rows() == reorderedB.rows())
    {
        for (size_t i = 0; i < B.rows(); ++i)
        {
            for (size_t j = 0; j < B.cols(); ++j)
            {
                if (B(i,j) == 0)
                {
                    BOOST_CHECK_SMALL(reorderedB(i,j), 1e-6);
                }
                else
                {
                    BOOST_CHECK_CLOSE_FRACTION(B(i,j), reorderedB(i,j), 1e-6);
                }
            }      
        }
    }
    else
    {
        std::cout << "Sizes don't match" << std::endl;
    }
}


BOOST_AUTO_TEST_CASE( TestFactorization ) {
    CompressorFixture compressor(quadrature);

    MatrixXd Q = compressor.getQ();
    MatrixXd R_11 = compressor.getR_11();

    MatrixXd B = compressor.getB();
    MatrixXd calculatedB = Q * R_11;
    std::vector<size_t> selectedK = compressor.getSelectedK();

    MatrixXd selectedB(Q.rows(), Q.cols());

    for (size_t i = 0; i < selectedK.size(); ++i)
    {
        selectedB.col(i) = B.col(selectedK[i]);
    }
}



BOOST_AUTO_TEST_SUITE_END()