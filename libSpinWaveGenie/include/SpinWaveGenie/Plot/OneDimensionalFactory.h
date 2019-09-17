//
//  OneDimensionalFactory.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/10/14.
//
//

#ifndef __OneDimensionalFactory__
#define __OneDimensionalFactory__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT OneDimensionalFactory
{
public:
  std::unique_ptr<OneDimensionalShapes> getGaussian(double fwhm, double tol);
  std::unique_ptr<OneDimensionalShapes> getGaussian(const std::array<double, 4> &fwhm, double tol);
  std::unique_ptr<OneDimensionalShapes> getLorentzian(double fwhm, double tol);
  std::unique_ptr<OneDimensionalShapes> getPseudoVoigt(double eta, double fwhm, double tol);

private:
};
} // namespace SpinWaveGenie
#endif /* defined(__OneDimensionalFactory__) */
