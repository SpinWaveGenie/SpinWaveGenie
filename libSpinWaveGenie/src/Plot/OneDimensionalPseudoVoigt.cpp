//
//  OneDimensionalPseudoVoigt.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/6/14.
//
//

#include <cmath>
#include "SpinWaveGenie/Plot/OneDimensionalPseudoVoigt.h"
#include "SpinWaveGenie/Plot/OneDimensionalFactory.h"

using namespace std;

namespace SpinWaveGenie
{

OneDimensionalPseudoVoigt::OneDimensionalPseudoVoigt() : eta(0.0), tolerance(0.001)
{
  Gaussian.reset(new OneDimensionalGaussian);
  Lorentzian.reset(new OneDimensionalLorentzian);
}

void OneDimensionalPseudoVoigt::setEta(double InEta) { eta = InEta; }

void OneDimensionalPseudoVoigt::setFWHM(double InFWHM)
{
  OneDimensionalFactory factory;
  Gaussian = factory.getGaussian(InFWHM, tolerance);
  Lorentzian = factory.getLorentzian(InFWHM, tolerance);
}

void OneDimensionalPseudoVoigt::setTolerance(double InTolerance)
{
  tolerance = InTolerance;
  Gaussian->setTolerance(InTolerance);
  Lorentzian->setTolerance(InTolerance);
}

double OneDimensionalPseudoVoigt::getMinimumEnergy() { return -1.0 * getMaximumEnergy(); }

double OneDimensionalPseudoVoigt::getMaximumEnergy() { return Lorentzian->getMaximumEnergy(); }

double OneDimensionalPseudoVoigt::getFunction(double frequency, double energy)
{
  return eta * Lorentzian->getFunction(frequency, energy) + (1.0 - eta) * Gaussian->getFunction(frequency, energy);
}

unique_ptr<OneDimensionalShapes> OneDimensionalPseudoVoigt::clone()
{
  return unique_ptr<OneDimensionalShapes>(new OneDimensionalPseudoVoigt(*this));
}
}
