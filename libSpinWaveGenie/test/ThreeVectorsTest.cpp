#define BOOST_TEST_MODULE EnergiesTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

using std::complex;
using namespace SpinWaveGenie;


BOOST_AUTO_TEST_CASE( InsertTest )
{
    ThreeVectors<double> doubleTest;
    ThreeVectors<complex<double>> complexTest;
    complex<double> complexNumber(1.0,1.0);
    doubleTest.insert(1.0,1.0,1.0);
    complexTest.insert(complexNumber,complexNumber,complexNumber);
    
    BOOST_CHECK(doubleTest.size() == 1);
    BOOST_CHECK(complexTest.size() == 1);

    const auto &point = *doubleTest.cbegin();
    BOOST_CHECK_CLOSE(point[0], 1.0, 1.0e-5);
    BOOST_CHECK_CLOSE(point[1], 1.0, 1.0e-5);
    BOOST_CHECK_CLOSE(point[2], 1.0, 1.0e-5);

    const auto &point2 = *complexTest.begin();
    BOOST_CHECK_SMALL(std::abs(point2[0] - complexNumber), 1.0e-5);
    BOOST_CHECK_SMALL(std::abs(point2[1] - complexNumber), 1.0e-5);
    BOOST_CHECK_SMALL(std::abs(point2[2] - complexNumber), 1.0e-5);
}

BOOST_AUTO_TEST_CASE( ClearTest )
{
    ThreeVectors<double> doubleTest;
    ThreeVectors<complex<double>> complexTest;
    complex<double> complexNumber(1.0,1.0);
    doubleTest.insert(1.0,1.0,1.0);
    complexTest.insert(complexNumber,complexNumber,complexNumber);
    BOOST_CHECK(doubleTest.size() == 1);
    BOOST_CHECK(complexTest.size() == 1);

    doubleTest.clear();
    complexTest.clear();

    BOOST_CHECK(doubleTest.size() == 0);
    BOOST_CHECK(complexTest.size() == 0);
}

BOOST_AUTO_TEST_CASE( IteratorTest )
{
    ThreeVectors<double> doubleTest;
    ThreeVectors<complex<double>> complexTest;
    complex<double> complexNumber(1.0,1.0);

    doubleTest.insert(1.0,1.0,1.0);
    doubleTest.insert(2.0,2.0,2.0);
    complexTest.insert(complexNumber,complexNumber,complexNumber);
    complexTest.insert(complexNumber*2.0,complexNumber*2.0,complexNumber*2.0);

    BOOST_CHECK(doubleTest.size() == 2);
    BOOST_CHECK(complexTest.size() == 2);

    int index = 0;
    for (const auto &point : doubleTest)
    {
        if (index == 0)
          BOOST_CHECK_CLOSE(point[0] + point[1] + point[2], 3.0, 1.0e-5);
        if (index == 1)
          BOOST_CHECK_CLOSE(point[0] + point[1] + point[2], 6.0, 1.0e-5);
        index++;
    }

    index = 0;
    for (const auto &point : complexTest)
    {
        if (index == 0)
          BOOST_CHECK_SMALL(std::abs(point[0] + point[1] + point[2] - complex<double>(3.0, 3.0)), 1.0e-5);
        if (index == 1)
          BOOST_CHECK_SMALL(std::abs(point[0] + point[1] + point[2] - complex<double>(6.0, 6.0)), 1.0e-5);
        index++;
    }

}
