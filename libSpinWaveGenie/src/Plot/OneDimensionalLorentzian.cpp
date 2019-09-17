//
//  OneDimensionalLorentzian.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//
#include "SpinWaveGenie/Plot/OneDimensionalLorentzian.h"

#include "boost/math/special_functions/pow.hpp"

#include <cmath>

namespace SpinWaveGenie
{

void OneDimensionalLorentzian::setFWHM(double InFWHM) { FWHM = InFWHM; }

void OneDimensionalLorentzian::setTolerance(double InTolerance) { Tolerance = InTolerance; }

double OneDimensionalLorentzian::getMinimumEnergy() { return -1.0 * getMaximumEnergy(); }

double OneDimensionalLorentzian::getMaximumEnergy()
{
  return (0.5 * FWHM) * sqrt(2.0 / (M_PI * FWHM * Tolerance) - 1.0);
}

double OneDimensionalLorentzian::getFunction(double energy)
{
  return FWHM / (2.0 * M_PI * (boost::math::pow<2>(m_frequency - energy) + boost::math::pow<2>(0.5 * FWHM)));
}

void OneDimensionalLorentzian::setFrequency(double frequency) { m_frequency = frequency; }

std::unique_ptr<OneDimensionalShapes> OneDimensionalLorentzian::clone() const
{
  return std::make_unique<OneDimensionalLorentzian>(*this);
}
} // namespace SpinWaveGenie
