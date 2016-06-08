#define BOOST_TEST_MODULE InteractionFactoryTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "SpinWaveGenie/SpinWaveGenie.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE(GetExchange)
{
  InteractionFactory factory;
  auto exchange = factory.getExchange("J", 1.0, "SA", "SB", 2.0, 2.1);
  BOOST_CHECK_EQUAL(exchange->getName(), "J");
  auto sl = exchange->sublattices();
  BOOST_CHECK_EQUAL(sl[0], "SA");
  BOOST_CHECK_EQUAL(sl[1], "SB");
}

BOOST_AUTO_TEST_CASE(GetDzyaloshinskiiMoriya)
{
  InteractionFactory factory;
  auto dm = factory.getDzyaloshinskiiMoriya("D", 1.0, {1.0, 0.0, 0.0}, "SA", "SB", 2.0, 2.1);
  BOOST_CHECK_EQUAL(dm->getName(), "D");
  auto sl = dm->sublattices();
  BOOST_CHECK_EQUAL(sl[0], "SA");
  BOOST_CHECK_EQUAL(sl[1], "SB");
}

BOOST_AUTO_TEST_CASE(getAnisotropy)
{
  InteractionFactory factory;
  auto anisotropy = factory.getAnisotropy("K", 1.0, {0.0, 0.0, 1.0}, "SA");
  BOOST_CHECK_EQUAL(anisotropy->getName(), "K");
  auto sl = anisotropy->sublattices();
  BOOST_CHECK_EQUAL(sl[0], "SA");
  BOOST_CHECK_EQUAL(sl[1], "SA");
}

BOOST_AUTO_TEST_CASE(GetMagneticField)
{
  InteractionFactory factory;
  auto magnetic_field = factory.getMagneticField("B", 1.0, {0.0, 0.0, 1.0}, "SA");
  BOOST_CHECK_EQUAL(magnetic_field->getName(), "B");
  auto sl = magnetic_field->sublattices();
  BOOST_CHECK_EQUAL(sl[0], "SA");
  BOOST_CHECK_EQUAL(sl[1], "SA");
}