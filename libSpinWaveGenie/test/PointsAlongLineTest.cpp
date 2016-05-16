#define BOOST_TEST_MODULE PointsAlongLineTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Containers/PointsAlongLine.h"

using namespace SpinWaveGenie;

ThreeVectors<double> getThreeVectors(int number)
{
    PointsAlongLine test;
    test.setFirstPoint(0.0,0.0,0.0);
    test.setFinalPoint(1.0,1.0,1.0);
    test.setNumberPoints(number);
    return test.getPoints();
}

BOOST_AUTO_TEST_CASE( ElevenPoints )
{
    double eps = 1.0e-5;
    ThreeVectors<double> test = getThreeVectors(11);
    BOOST_CHECK(test.size() == 11);
    
    double value = 0.0;
    for (const auto &point : test)
    {
      BOOST_CHECK_CLOSE(value, point[0], eps);
      BOOST_CHECK_CLOSE(value, point[1], eps);
      BOOST_CHECK_CLOSE(value, point[2], eps);
      value = value + 0.1;
    }
}

BOOST_AUTO_TEST_CASE( OnePoint )
{
    double eps = 1.0e-5;
    ThreeVectors<double> test = getThreeVectors(1);
    BOOST_CHECK(test.size() == 1);
    
    double value = 0.0;
    for (const auto &point : test)
    {
      BOOST_CHECK_CLOSE(value, point[0], eps);
      BOOST_CHECK_CLOSE(value, point[1], eps);
      BOOST_CHECK_CLOSE(value, point[2], eps);
    }
}
