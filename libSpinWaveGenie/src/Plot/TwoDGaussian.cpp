//
//  TwoDGaussian.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/23/14.
//
//
#define NOMINMAX
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

double TwoDGaussian::getFunction(double frequency, double energy)
{
  return exp(-1.0 * (c * pow(frequency - energy, 2) + 2.0 * b * (frequency - energy) * u + a * pow(u, 2)));
}

unique_ptr<OneDimensionalShapes> TwoDGaussian::clone()
{
  return unique_ptr<OneDimensionalShapes>(new TwoDGaussian(*this));
}
}