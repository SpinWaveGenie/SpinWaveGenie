//
//  OneDimensionalPseudoVoigt.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/6/14.
//
//

#include "SpinWaveGenie/Plot/OneDimensionalPseudoVoigt.h"
#include "SpinWaveGenie/Plot/OneDimensionalFactory.h"
#include <cmath>

using namespace std;

namespace SpinWaveGenie
{

OneDimensionalPseudoVoigt::OneDimensionalPseudoVoigt() = default;

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

void OneDimensionalPseudoVoigt::setFrequency(double frequency)
{
  Lorentzian->setFrequency(frequency);
  Gaussian->setFrequency(frequency);
}

double OneDimensionalPseudoVoigt::getFunction(double energy)
{
  return eta * Lorentzian->getFunction(energy) + (1.0 - eta) * Gaussian->getFunction(energy);
}

unique_ptr<OneDimensionalShapes> OneDimensionalPseudoVoigt::clone() const
{
  return std::make_unique<OneDimensionalPseudoVoigt>(*this);
}
} // namespace SpinWaveGenie
