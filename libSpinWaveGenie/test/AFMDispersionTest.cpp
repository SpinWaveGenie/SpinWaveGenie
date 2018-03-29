#define BOOST_TEST_MODULE AFMDispersionTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include "SpinWaveGenie/SpinWaveGenie.h"

using namespace std;
using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE(DispersionTest)
{
  double J = 1.0;
  double D = 0.1;
  double S = 0.5;
  Cell cell;
  cell.setBasisVectors(1.0, 10.0, 10.0, 90.0, 90.0, 90.0);

  Sublattice spin0;
  string name0 = "Spin0";
  spin0.setName(name0);
  spin0.setType("NONE");
  spin0.setMoment(S, 0.0, 0.0);
  cell.addSublattice(spin0);
  cell.addAtom(name0, 0.0, 0.0, 0.0);

  Sublattice spin1;
  string name1 = "Spin1";
  spin1.setName(name1);
  spin1.setType("NONE");
  spin1.setMoment(S, M_PI, 0.0);
  cell.addSublattice(spin1);
  cell.addAtom(name1, 0.5, 0.0, 0.0);

  SpinWaveBuilder builder(cell);

  InteractionFactory interactions;

  Eigen::Vector3d xhat(1.0, 0.0, 0.0);
  builder.addInteraction(interactions.getExchange("J", -1.0 * J, name0, name1, 0.4, 0.6));
  builder.addInteraction(interactions.getAnisotropy("D", D, xhat, name0));
  builder.addInteraction(interactions.getAnisotropy("D", D, xhat, name1));

  SpinWave genie = builder.createElement();

  PointsAlongLine Line;
  Line.setFirstPoint(0.0, 0.0, 0.0);
  Line.setFinalPoint(3.0, 0.0, 0.0);
  Line.setNumberPoints(64);
  ThreeVectors<double> kPoints = Line.getPoints();

  for (const auto &kpt : kPoints)
  {
    double k = kpt[0];
    genie.createMatrix(k, 0.0, 0.0);
    genie.calculate();
    Results pts = genie.getPoints();
    pts.sort();

    // analytical solution for frequency and intensity
    Results soln;
    Point tmp;
    tmp.frequency = sqrt(4.0 * J * S * S * (1.0 - cos(M_PI * k)) * (J * (1.0 + cos(M_PI * k)) + D));
    tmp.intensity = 0.25 * S * sqrt(J * (1.0 - cos(M_PI * k)) / (J * (1 + cos(M_PI * k)) + D));
    soln.insert(tmp);

    tmp.frequency = sqrt(4.0 * J * S * S * (1.0 + cos(M_PI * k)) * (J * (1.0 - cos(M_PI * k)) + D));
    double denominator = (J * (1.0 + cos(M_PI * k)));
    if (std::abs(denominator) < std::numeric_limits<double>::epsilon()) {
      tmp.intensity = std::numeric_limits<double>::quiet_NaN();
    } else {
      tmp.intensity = 0.25 * S * sqrt((J * (1.0 - cos(M_PI * k)) + D) / denominator);
}
    soln.insert(tmp);
    soln.sort();

    // intensity return nan when the frequency is 0.0
    if (std::abs(pts[0].frequency) > 1.0e-5)
    {
      BOOST_CHECK_CLOSE(pts[0].frequency, soln[0].frequency, 1.0e-5);
      BOOST_CHECK_CLOSE(pts[1].frequency, soln[1].frequency, 1.0e-5);
      if (cos(k * M_PI) > 0.0) {
        BOOST_CHECK_CLOSE(pts[0].intensity + pts[1].intensity, soln[1].intensity, 1.0e-5);
      } else {
        BOOST_CHECK_CLOSE(pts[0].intensity + pts[1].intensity, soln[0].intensity, 1.0e-5);
}
    }
  }
}
