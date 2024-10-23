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


double lowerBound = 0;
double upperBound = 2; 
std::vector<double> param1 = {5, 4};
std::vector<double> param2 = {6, 3};
FunctionHandler<std::vector<double>, std::vector<double>> handlerPiecewiseSmooth(testFunction2ParamPC, param1, param2);

QuadratureRule<FunctionHandler<std::vector<double>, std::vector<double>>> quadrature(lowerBound, upperBound, handlerPiecewiseSmooth);

template<typename T>
struct CompressorFixture: public Compressor<T>
{
    CompressorFixture(const T& quadratureInput) : Compressor<T>(quadratureInput) {}

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
    quadrature.calculateQuadratureNodes();
    CompressorFixture compressor(quadrature);

    matrix<double> A = compressor.getA();

    BOOST_CHECK_EQUAL(A.size1(), 60);
    BOOST_CHECK_EQUAL(A.size2(), 4);

    matrix_column<matrix<double>> col1(A, 0);

    for (std::size_t i = 0; i < col1.size(); ++i) {
        std::cout << col1(i) << std::endl;
    }
}

BOOST_AUTO_TEST_SUITE_END()