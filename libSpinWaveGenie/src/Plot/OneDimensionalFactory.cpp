//
//  OneDimensionalFactory.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/10/14.
//
//

#include "SpinWaveGenie/Plot/OneDimensionalFactory.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/OneDimensionalLorentzian.h"
#include "SpinWaveGenie/Plot/OneDimensionalPseudoVoigt.h"

namespace SpinWaveGenie
{

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getGaussian(double fwhm, double tol)
{
  std::unique_ptr<OneDimensionalGaussian> resinfo(new OneDimensionalGaussian);
  resinfo->setFWHM(fwhm);
  resinfo->setTolerance(tol);
  return std::move(resinfo);
}

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getLorentzian(double fwhm, double tol)
{
  std::unique_ptr<OneDimensionalLorentzian> resinfo(new OneDimensionalLorentzian);
  resinfo->setFWHM(fwhm);
  resinfo->setTolerance(tol);
  return std::move(resinfo);
}

std::unique_ptr<OneDimensionalShapes> OneDimensionalFactory::getPseudoVoigt(double eta, double fwhm, double tol)
{
  std::unique_ptr<OneDimensionalPseudoVoigt> resinfo(new OneDimensionalPseudoVoigt);
  resinfo->setEta(eta);
  resinfo->setFWHM(fwhm);
  resinfo->setTolerance(tol);
  return std::move(resinfo);
}
}