//
//  File.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/30/14.
//
//

#ifndef __OneDimensionalGaussian__
#define __OneDimensionalGaussian__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT OneDimensionalGaussian : public OneDimensionalShapes
{
public:
  OneDimensionalGaussian(double FWHM = 1.0, double Tolerance = 0.01);
  void setFWHM(double InFWHM);
  void setTolerance(double InTolerance) override;
  double getMinimumEnergy() override;
  double getMaximumEnergy() override;
  double getFunction(double frequency, double energy) override;
  std::unique_ptr<OneDimensionalShapes> clone() override;

private:
  void update();
  double m_FWHM, m_Tolerance;
  double m_Diff, m_Factor, m_ma;
};
}
#endif /* defined(__OneDimensionalGaussian__) */
