//
//  IntegrateThetaPhi.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#include <vector>
#include "SpinWaveGenie/Plot/IntegrateThetaPhi.h"
#include "SpinWaveGenie/Plot/AdaptiveSimpson.h"

using namespace std;

namespace SpinWaveGenie
{



std::unique_ptr<SpinWavePlot> IntegrateThetaPhi::clone()
{
  return unique_ptr<SpinWavePlot>(new IntegrateThetaPhi(*this));
}

const Cell &IntegrateThetaPhi::getCell() const { return resolutionFunction->getCell(); }

const Energies &IntegrateThetaPhi::getEnergies() { return resolutionFunction->getEnergies(); }

void IntegrateThetaPhi::setEnergies(Energies energiesIn) { resolutionFunction->setEnergies(energiesIn); }

std::vector<double> IntegrateThetaPhi::calculateIntegrand(std::deque<double> &x)
{
  Vector3 tmp, k;
  double theta = x[0];
  double phi = x[1];

  // cout << "r= " << r << endl;
  // cout << "theta= " << theta << endl;
  // cout << "phi= " << phi << endl;

  tmp[0] = r * sin(theta) * cos(phi);
  tmp[1] = r * sin(theta) * sin(phi);
  tmp[2] = r * cos(theta);

  Matrix3 basisVectors = resolutionFunction->getCell().getBasisVectors();

  k = tmp.transpose() * basisVectors / (2.0 * M_PI);
  // cout << tmp.norm() << endl;
  // cout << k.transpose() << endl;

  vector<double> val = resolutionFunction->getCut(k[0], k[1], k[2]);

  // for (auto it = values.begin(); it!=values.end(); it++)
  //{
  //    cout << (*it) << endl;
  //}

  // cout << MinimumEnergy << " " << MaximumEnergy << " " << EnergyPoints << endl;
  double factor = sin(theta) / (4.0 * M_PI);
  std::transform(val.begin(), val.end(), val.begin(),
                 std::bind(std::multiplies<double>(), factor, std::placeholders::_1));
  return val;
}

std::vector<double> IntegrateThetaPhi::getCut(double kx, double ky, double kz)
{
  std::vector<double> xmin = {0.0, 0.0};
  std::vector<double> xmax = {M_PI, 2.0 * M_PI};
  r = std::abs(kz);

  // cout << "dispAng = " << r << endl;

  auto funct = std::bind<std::vector<double>>(&IntegrateThetaPhi::calculateIntegrand, this, std::placeholders::_1);
  // GaussKronrod test;
  AdaptiveSimpson test;
  // AdaptiveBoole test;
  test.setFunction(funct);
  test.setInterval(xmin, xmax);
  test.setPrecision(tolerance);
  test.setMaximumDivisions(maximumEvaluations);
  return test.integrate();
}
}
