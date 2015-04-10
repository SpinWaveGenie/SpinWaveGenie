#define BOOST_TEST_MODULE FMDispersionTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include "SpinWaveGenie/SpinWaveGenie.h"

using namespace std;
using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE(DispersionTest)
{
  Cell cell;
  cell.setBasisVectors(1.0, 10.0, 10.0, 90.0, 90.0, 90.0);

  double J = 1.0;
  double S = 1.0;

  Sublattice spin0;
  string name0 = "Spin0";
  spin0.setName(name0);
  spin0.setType("NONE");
  spin0.setMoment(S, 0.0, 0.0);
  cell.addSublattice(spin0);
  cell.addAtom(name0, 0.0, 0.0, 0.0);

  SpinWaveBuilder builder(cell);
  InteractionFactory interactions;
  builder.addInteraction(interactions.getExchange("J", J, name0, name0, 0.9, 1.1));
  SpinWave genie = builder.createElement();

  PointsAlongLine line;
  line.setFirstPoint(0.0, 0.0, 0.0);
  line.setFinalPoint(3.0, 0.0, 0.0);
  line.setNumberPoints(61);
  ThreeVectors<double> kPoints = line.getPoints();

  for (auto it = kPoints.begin(); it != kPoints.end(); ++it)
  {
    double k = it->get<0>();
    genie.createMatrix(k, 0.0, 0.0);
    genie.calculate();
    Results pt = genie.getPoints();

    //analytical solution for frequency and intensity
    double frequency = 2.0 * J * S * (1.0 - cos(2.0 * M_PI * k));
    double intensity = S / 4.0;

    BOOST_CHECK_CLOSE(pt[0].frequency, frequency, 1.0e-5);

    // intensity return nan when the frequency is 0.0
    if (std::abs(pt[0].frequency) > 1.0e-5)
    {
      BOOST_CHECK_CLOSE(pt[0].intensity, intensity, 1.0e-5);
    }
  }
}
