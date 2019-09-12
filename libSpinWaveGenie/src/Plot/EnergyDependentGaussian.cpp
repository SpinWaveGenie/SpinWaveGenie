#include "SpinWaveGenie/Plot/EnergyDependentGaussian.h"

#include "boost/math/special_functions/pow.hpp"

#include <cmath>
#include <limits>

namespace SpinWaveGenie
{

void EnergyDependentGaussian::update()
{
  constexpr double factor = 0.5 / std::sqrt(2. * std::log(2.));
  for (auto &value : m_sigma)
  {
    value *= factor;
  }
}

EnergyDependentGaussian::EnergyDependentGaussian(const std::array<double, 4> &FWHM, double Tolerance)
    : m_sigma(FWHM), m_Tolerance(Tolerance)
{
  this->update();
}

void EnergyDependentGaussian::setFWHM(const std::array<double, 4> &InFWHM)
{
  m_sigma = InFWHM;
  this->update();
}

void EnergyDependentGaussian::setTolerance(double InTolerance)
{
  m_Tolerance = InTolerance;
  this->update();
}

double EnergyDependentGaussian::getMinimumEnergy() { return std::numeric_limits<double>::lowest(); }

double EnergyDependentGaussian::getMaximumEnergy() { return std::numeric_limits<double>::max(); }

double EnergyDependentGaussian::getFunction(double frequency, double energy)
{
  double sigma = m_sigma[0] + m_sigma[1] * frequency + m_sigma[2] * boost::math::pow<2>(frequency) +
                 m_sigma[3] * boost::math::pow<3>(frequency);
  std::cout << "sigma: " << sigma << '\n';
  return 1.0 / (sigma * std::sqrt(2.0 * M_PI)) * std::exp(-1.0 / (2.0 * boost::math::pow<2>(sigma)) *
                boost::math::pow<2>(frequency - energy));
}

std::unique_ptr<OneDimensionalShapes> EnergyDependentGaussian::clone() const
{
  return std::make_unique<EnergyDependentGaussian>(*this);
}
} // namespace SpinWaveGenie
