#include "SpinWaveGenie/Plot/EnergyDependentGaussian.h"

#include "boost/math/special_functions/pow.hpp"

#include <cmath>
#include <limits>

namespace SpinWaveGenie
{

void EnergyDependentGaussian::update()
{
  double factor = 2.0 * sqrt(log(2.0)) / (m_FWHM * sqrt(M_PI));
  double constant = log(m_Tolerance / factor);
  m_ma = -4.0 * log(2.0) / pow(m_FWHM, 2);
  m_Diff = sqrt(constant / m_ma);
  m_Factor = 2.0 * sqrt(log(2.0)) / (m_FWHM * sqrt(M_PI));
}

EnergyDependentGaussian::EnergyDependentGaussian(const std::array<double, 4> &FWHM, double Tolerance) : m_expansion(FWHM), m_Tolerance(Tolerance)
{
}

void EnergyDependentGaussian::setFWHM(const std::array<double, 4> &InFWHM)
{
  m_expansion = InFWHM;
}

void EnergyDependentGaussian::setFWHM(double frequency) {
  m_FWHM = m_expansion[0] + m_expansion[1] * frequency + m_expansion[2] * boost::math::pow<2>(frequency) + m_expansion[3] * boost::math::pow<3>(frequency);
  this->update();
}

void EnergyDependentGaussian::setTolerance(double InTolerance)
{
  m_Tolerance = InTolerance;
}

double EnergyDependentGaussian::getMinimumEnergy() { return -1.0 * m_Diff; }

double EnergyDependentGaussian::getMaximumEnergy()
{
  return m_Diff;
}

double EnergyDependentGaussian::getFunction(double frequency, double energy)
{
  this->setFWHM(frequency);
  return m_Factor * exp(m_ma * pow(frequency - energy, 2));
}

std::unique_ptr<OneDimensionalShapes> EnergyDependentGaussian::clone() const
{
  return std::make_unique<EnergyDependentGaussian>(*this);
}
} // namespace SpinWaveGenie
