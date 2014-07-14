#define BOOST_TEST_MODULE HKLDirectionsTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Containers/HKLDirections.h"


using namespace SpinWaveGenie;


bool areEqual(Axis first, Axis second)
{
    if (first.v0 == second.v0 &&
        first.v1 == second.v1 &&
        first.v2 == second.v2 &&
        first.delta == second.delta)
        return true;
    else
        return false;
}

HKLDirections getDirections()
{
    HKLDirections directions;
    directions.addDirection(0, 0.05);
    directions.addDirection(0.5,-1.0,0.0,0.1);
    directions.addDirection(2, 0.05);
    return directions;
}


BOOST_AUTO_TEST_CASE( iteratorAccess )
{
    Axis test1,test2,test3;
    test1.v0 = 1.0;
    test1.v1 = 0.0;
    test1.v2 = 0.0;
    test1.delta = 0.05;
    
    test2.v0 =  0.5;
    test2.v1 = -1.0;
    test2.v2 =  0.0;
    test2.delta = 0.1;
    
    test3.v0 =  0.0;
    test3.v1 =  0.0;
    test3.v2 =  1.0;
    test3.delta = 0.05;
    
    HKLDirections directions = getDirections();
    
    BOOST_CHECK(directions.size() == 3);
    
    auto it = directions.begin();
    
    BOOST_CHECK(areEqual(test1,*it));
    it++;
    BOOST_CHECK(areEqual(test2,*it));
    it++;
    BOOST_CHECK(areEqual(test3,*it));
}

BOOST_AUTO_TEST_CASE( subscriptAccess )
{
    Axis test1,test2,test3;
    test1.v0 = 1.0;
    test1.v1 = 0.0;
    test1.v2 = 0.0;
    test1.delta = 0.05;
    
    test2.v0 =  0.5;
    test2.v1 = -1.0;
    test2.v2 =  0.0;
    test2.delta = 0.1;
    
    test3.v0 =  0.0;
    test3.v1 =  0.0;
    test3.v2 =  1.0;
    test3.delta = 0.05;
    
    HKLDirections directions = getDirections();
    
    BOOST_CHECK(directions.size() == 3);
    
    BOOST_CHECK(areEqual(test1,directions[0]));
    BOOST_CHECK(areEqual(test2,directions[1]));
    BOOST_CHECK(areEqual(test3,directions[2]));
}

