//
//  IntegrateEnergy.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#include "SpinWaveGenie/Plot/IntegrateEnergy.h"
#include "SpinWaveGenie/Plot/AdaptiveSimpson.h"

using namespace std;

namespace SpinWaveGenie
{

IntegrateEnergy::IntegrateEnergy(const IntegrateEnergy &other) : kx(0.0),ky(0.0),kz(0.0)
{
  resolutionFunction = move(other.resolutionFunction->clone());
  this->centeredEnergies = other.centeredEnergies;
  this->delta = other.delta;
  this->tolerance = other.tolerance;
  this->maximumEvaluations = other.maximumEvaluations;
}

IntegrateEnergy::IntegrateEnergy(unique_ptr<SpinWavePlot> resFunction, Energies energies, double delta, double tol,
                                 int maxEvals)
    : resolutionFunction(std::move(resFunction)), maximumEvaluations(maxEvals), tolerance(tol), delta(delta), kx(0.0),
      ky(0.0), kz(0.0), centeredEnergies(energies)
{
}

std::vector<double> IntegrateEnergy::calculateIntegrand(std::deque<double> &x)
{
  assert(x.size() == 1);

  Energies newEnergies;
  for (auto value = centeredEnergies.begin(); value != centeredEnergies.end(); value++)
  {
    newEnergies.insert(*value + x[0]);
  }

  resolutionFunction->setEnergies(newEnergies);

  return resolutionFunction->getCut(kx, ky, kz);
}

std::vector<double> IntegrateEnergy::getCut(double kxIn, double kyIn, double kzIn)
{
  kx = kxIn;
  ky = kyIn;
  kz = kzIn;

  auto funct = std::bind<std::vector<double>>(&IntegrateEnergy::calculateIntegrand, this, std::placeholders::_1);
  AdaptiveSimpson test;
  test.setFunction(funct);
  std::vector<double> xmin = {-1.0 * delta};
  std::vector<double> xmax = {delta};
  test.setInterval(xmin, xmax);
  test.setPrecision(tolerance);
  test.setMaximumDivisions(maximumEvaluations);
  return test.integrate();
}

const Cell &IntegrateEnergy::getCell() const { return resolutionFunction->getCell(); }

const Energies &IntegrateEnergy::getEnergies() { return centeredEnergies; }

void IntegrateEnergy::setEnergies(Energies energiesIn) { centeredEnergies = energiesIn; }

std::unique_ptr<SpinWavePlot> IntegrateEnergy::clone() { return unique_ptr<SpinWavePlot>(new IntegrateEnergy(*this)); }
}
