#define BOOST_TEST_MODULE DGNTestSuite  
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/algorithm/cxx11/is_sorted.hpp>
#include "ggq.hpp"
#include "damped_gauss_newton.hpp"


struct DGNFixture: public DGN 
{};


BOOST_AUTO_TEST_SUITE( DGNTestSuite )

BOOST_AUTO_TEST_SUITE_END()