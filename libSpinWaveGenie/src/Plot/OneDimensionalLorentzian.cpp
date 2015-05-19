//
//  OneDimensionalLorentzian.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//
#define _USE_MATH_DEFINES
#include <cmath>
#include "SpinWaveGenie/Plot/OneDimensionalLorentzian.h"

using namespace std;

namespace SpinWaveGenie
{

void OneDimensionalLorentzian::setFWHM(double InFWHM) { FWHM = InFWHM; }

void OneDimensionalLorentzian::setTolerance(double InTolerance) { Tolerance = InTolerance; }

double OneDimensionalLorentzian::getMinimumEnergy() { return -1.0 * getMaximumEnergy(); }

double OneDimensionalLorentzian::getMaximumEnergy()
{
  return (0.5 * FWHM) * sqrt(2.0 / (M_PI * FWHM * Tolerance) - 1.0);
}

double OneDimensionalLorentzian::getFunction(double frequency, double energy)
{
  return FWHM / (2.0 * M_PI * (pow(frequency - energy, 2) + pow(0.5 * FWHM, 2)));
}

unique_ptr<OneDimensionalShapes> OneDimensionalLorentzian::clone()
{
  return unique_ptr<OneDimensionalShapes>(new OneDimensionalLorentzian(*this));
}
}