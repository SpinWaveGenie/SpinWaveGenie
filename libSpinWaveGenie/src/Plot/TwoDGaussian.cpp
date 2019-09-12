//
//  TwoDGaussian.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/23/14.
//
//
#include <cmath>
#include <algorithm>
#include "SpinWaveGenie/Plot/TwoDGaussian.h"

using namespace std;

namespace SpinWaveGenie
{

void TwoDGaussian::setResolution(double aIn, double bIn, double cIn)
{
  a = aIn;
  b = bIn;
  c = cIn;
}

void TwoDGaussian::setU(double uIn) { u = uIn; }

void TwoDGaussian::setTolerance(double InTolerance) { Tolerance = InTolerance; }

double TwoDGaussian::getMinimumEnergy()
{
  double d = -1.0 * log(Tolerance);
  double firstSolution = (-b * u + sqrt((b * b - a * c) * u * u + c * d)) / c;
  double secondSolution = (a * u * u - d) / (c * firstSolution);
  return min(firstSolution, secondSolution);
}

double TwoDGaussian::getMaximumEnergy()
{
  double d = -1.0 * log(Tolerance);
  double firstSolution = (-b * u + sqrt((b * b - a * c) * u * u + c * d)) / c;
  double secondSolution = (a * u * u - d) / (c * firstSolution);
  return max(firstSolution, secondSolution);
}

void TwoDGaussian::setFrequency(double frequency)
{
    m_frequency = frequency;
}

double TwoDGaussian::getFunction(double energy)
{
  return exp(-1.0 * (c * pow(m_frequency - energy, 2) + 2.0 * b * (m_frequency - energy) * u + a * pow(u, 2)));
}

unique_ptr<OneDimensionalShapes> TwoDGaussian::clone() const { return std::make_unique<TwoDGaussian>(*this); }
}
