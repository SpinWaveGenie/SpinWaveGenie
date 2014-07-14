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

    auto it = doubleTest.begin();
    BOOST_CHECK_CLOSE(it->get<0>(),1.0,1.0e-5);
    BOOST_CHECK_CLOSE(it->get<1>(),1.0,1.0e-5);
    BOOST_CHECK_CLOSE(it->get<2>(),1.0,1.0e-5);

    auto it2 = complexTest.begin();
    BOOST_CHECK_SMALL(std::abs(it2->get<0>()-complexNumber),1.0e-5);
    BOOST_CHECK_SMALL(std::abs(it2->get<1>()-complexNumber),1.0e-5);
    BOOST_CHECK_SMALL(std::abs(it2->get<2>()-complexNumber),1.0e-5);
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
    for (auto it = doubleTest.begin(); it!=doubleTest.end(); ++it)
    {
        if (index == 0)
            BOOST_CHECK_CLOSE(it->get<0>()+it->get<1>()+it->get<2>(),3.0,1.0e-5);
        if (index == 1)
            BOOST_CHECK_CLOSE(it->get<0>()+it->get<1>()+it->get<2>(),6.0,1.0e-5);
        index++;
    }

    index = 0;
    for (auto it = complexTest.begin(); it!=complexTest.end(); ++it)
    {
        if (index == 0)
            BOOST_CHECK_SMALL(std::abs(it->get<0>()+it->get<1>()+it->get<2>()-complex<double>(3.0,3.0)),1.0e-5);
        if (index == 1)
            BOOST_CHECK_SMALL(std::abs(it->get<0>()+it->get<1>()+it->get<2>()-complex<double>(6.0,6.0)),1.0e-5);
        index++;
    }

}
