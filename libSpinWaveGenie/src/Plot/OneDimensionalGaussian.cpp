#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"

#include "boost/math/special_functions/pow.hpp"

#include <cmath>

using namespace std;

namespace SpinWaveGenie
{

void OneDimensionalGaussian::update()
{
  double factor = 2.0 * sqrt(log(2.0)) / (m_FWHM * sqrt(M_PI));
  double constant = log(m_Tolerance / factor);
  m_ma = -4.0 * log(2.0) / pow(m_FWHM, 2);
  m_Diff = sqrt(constant / m_ma);
  m_Factor = 2.0 * sqrt(log(2.0)) / (m_FWHM * sqrt(M_PI));
}

OneDimensionalGaussian::OneDimensionalGaussian(double FWHM, double Tolerance) : m_FWHM(FWHM), m_Tolerance(Tolerance)
{
  this->update();
}

void OneDimensionalGaussian::setFWHM(double InFWHM)
{
  m_FWHM = InFWHM;
  this->update();
}

void OneDimensionalGaussian::setTolerance(double InTolerance)
{
  m_Tolerance = InTolerance;
  this->update();
}

double OneDimensionalGaussian::getMinimumEnergy() { return -1.0 * m_Diff; }

double OneDimensionalGaussian::getMaximumEnergy()
{
  return m_Diff;
}

void OneDimensionalGaussian::setFrequency(double frequency)
{
    m_frequency = frequency;
}

double OneDimensionalGaussian::getFunction(double energy)
{
  return m_Factor * exp(m_ma * boost::math::pow<2>(m_frequency - energy));
}

unique_ptr<OneDimensionalShapes> OneDimensionalGaussian::clone() const
{
  return std::make_unique<OneDimensionalGaussian>(*this);
}
}
