//
//  File.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/30/14.
//
//

#ifndef __OneDimensionalGaussian__
#define __OneDimensionalGaussian__

#include <iostream>
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"

namespace SpinWaveGenie
{

class OneDimensionalGaussian : public OneDimensionalShapes
{
public:
  OneDimensionalGaussian(double FWHM = 1.0, double Tolerance = 0.01);
  void setFWHM(double InFWHM);
  void setTolerance(double InTolerance);
  double getMinimumEnergy();
  double getMaximumEnergy();
  double getFunction(double frequency, double energy);
  std::unique_ptr<OneDimensionalShapes> clone();
  ~OneDimensionalGaussian(){};

private:
  void update();
  double m_FWHM, m_Tolerance;
  double m_Diff, m_Factor, m_ma;
};
}
#endif /* defined(__OneDimensionalGaussian__) */
