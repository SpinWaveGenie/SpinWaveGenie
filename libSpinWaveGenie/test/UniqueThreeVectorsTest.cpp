#define BOOST_TEST_MODULE EnergiesTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE( InsertTest )
{
    UniqueThreeVectors<double> doubleTest;
    doubleTest.insert(1.0,1.0,1.0);
    doubleTest.insert(1.0,1.0,1.0);
    
    BOOST_CHECK(doubleTest.size() == 1);

    const auto &point = *doubleTest.begin();
    BOOST_CHECK_CLOSE(point[0], 1.0, 1.0e-5);
    BOOST_CHECK_CLOSE(point[1], 1.0, 1.0e-5);
    BOOST_CHECK_CLOSE(point[2], 1.0, 1.0e-5);
}

BOOST_AUTO_TEST_CASE( ClearTest )
{
    UniqueThreeVectors<double> doubleTest;
    doubleTest.insert(1.0,1.0,1.0);
    doubleTest.insert(1.0,1.0,1.0);
    BOOST_CHECK(doubleTest.size() == 1);

    doubleTest.clear();

    BOOST_CHECK(doubleTest.empty());
}

BOOST_AUTO_TEST_CASE( IteratorTest )
{
    UniqueThreeVectors<double> doubleTest;

    doubleTest.insert(1.0,1.0,1.0);
    doubleTest.insert(1.0,1.0,1.0);

    doubleTest.insert(2.0,2.0,2.0);
    doubleTest.insert(2.0,2.0,2.0);

    BOOST_CHECK(doubleTest.size() == 2);

    int index = 0;
    for (const auto &point : doubleTest)
    {
        if (index == 0)
          BOOST_CHECK_CLOSE(point[0] + point[1] + point[2], 3.0, 1.0e-5);
        if (index == 1)
          BOOST_CHECK_CLOSE(point[0] + point[1] + point[2], 6.0, 1.0e-5);
        index++;
    }
}
