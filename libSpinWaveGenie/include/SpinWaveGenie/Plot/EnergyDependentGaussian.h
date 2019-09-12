//
//  File.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/30/14.
//
//

#ifndef __EnergyDependentGaussian__
#define __EnergyDependentGaussian__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT EnergyDependentGaussian : public OneDimensionalShapes
{
public:
  EnergyDependentGaussian(const std::array<double, 4> &FWHM, double Tolerance = 0.01);
  void setFWHM(const std::array<double, 4> &FWHM);
  void setTolerance(double InTolerance) override;
  double getMinimumEnergy() override;
  double getMaximumEnergy() override;
  double getFunction(double frequency, double energy) override;
  std::unique_ptr<OneDimensionalShapes> clone() const override;

private:
  void update();
  void setFWHM(double frequency);
  std::array<double, 4> m_expansion;
  double m_FWHM, m_Tolerance;
  double m_Diff, m_Factor, m_ma;
};
} // namespace SpinWaveGenie
#endif /* defined(__EnergyDependentGaussian__) */
