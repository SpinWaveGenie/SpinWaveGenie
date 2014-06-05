#define BOOST_TEST_MODULE EnergiesTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "Containers/Energies.h"
#include <iostream>
#include <stdexcept>
#include <string>


BOOST_AUTO_TEST_CASE( ConstructorsTest )
{
    Energies energies(0.0,100.0,101);
    
    BOOST_CHECK(energies.size() == 101);
}

BOOST_AUTO_TEST_CASE( LowerUpperBoundTest )
{
    Energies energies(0.0,100.0,101);
    
    BOOST_CHECK(energies.getLowerBound(35.0) ==35);
    BOOST_CHECK(energies.getUpperBound(85.0) ==86);
    
    BOOST_CHECK(energies.getLowerBound(-10.0) == 0);
    BOOST_CHECK(energies.getUpperBound(-10.0) == 0);
    
    BOOST_CHECK(energies.getLowerBound(110.0) == 101);
    BOOST_CHECK(energies.getUpperBound(110.0) == 101);
}

BOOST_AUTO_TEST_CASE( subscriptOperatorTest )
{
    Energies energies(0.0,100.0,101);
    
    BOOST_CHECK_CLOSE(energies[35],35.0,1.0e-5);
    BOOST_CHECK_CLOSE(energies[85],85.0,1.0e-5);
    
}

BOOST_AUTO_TEST_CASE( clearTest )
{
    Energies energies(0.0,100.0,101);
    BOOST_CHECK(energies.size() == 101);
    energies.clear();
    BOOST_CHECK(energies.size() == 0);
}


BOOST_AUTO_TEST_CASE( iteratorTest )
{

    Energies energies(0.0,10.0,11);
    
    double value = 0.0;
    for (auto it = energies.begin();it!=energies.end();it++)
    {
        BOOST_CHECK_CLOSE(*it,value,1.0e-5);
        value += 1.0;
    }
}







