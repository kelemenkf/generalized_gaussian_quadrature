#define BOOST_TEST_MODULE CompressorTestSuite
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "compressor.hpp"


template<typename T>
struct CompressorFixture: public Compressor<T>
{
    CommpressorFixture() : Compressor<T>() {}

    ~CompressorFixture() {}
}; 


BOOST_AUTO_TEST_SUITE( CompressorTestSuite )


BOOST_AUTO_TEST_CASE( TestCompressorConstructor ) {

}

BOOST_AUTO_TEST_CASE( TestConstructA ) {

}

BOOST_AUTO_TEST_SUITE_END()