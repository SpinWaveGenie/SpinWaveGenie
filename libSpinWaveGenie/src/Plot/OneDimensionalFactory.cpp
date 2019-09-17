//
//  OneDimensionalFactory.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/10/14.
//
//

#include "SpinWaveGenie/Plot/OneDimensionalFactory.h"
#include "SpinWaveGenie/Plot/EnergyDependentGaussian.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/OneDimensionalLorentzian.h"
#include "SpinWaveGenie/Plot/OneDimensionalPseudoVoigt.h"

namespace SpinWaveGenie
{

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getGaussian(double fwhm, double tol)
{
  return std::make_unique<OneDimensionalGaussian>(fwhm, tol);
}

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getGaussian(const std::array<double,4> &fwhm, double tol)
{
  return std::make_unique<EnergyDependentGaussian>(fwhm, tol);
}

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getLorentzian(double fwhm, double tol)
{
  auto resinfo = std::make_unique<OneDimensionalLorentzian>();
  resinfo->setFWHM(fwhm);
  resinfo->setTolerance(tol);
  return resinfo;
}

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getPseudoVoigt(double eta, double fwhm, double tol)
{
  auto resinfo = std::make_unique<OneDimensionalPseudoVoigt>();
  resinfo->setEta(eta);
  resinfo->setFWHM(fwhm);
  resinfo->setTolerance(tol);
  return resinfo;
}
} // namespace SpinWaveGenie
