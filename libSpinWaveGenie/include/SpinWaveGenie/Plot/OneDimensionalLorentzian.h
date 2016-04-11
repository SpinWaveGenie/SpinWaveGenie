//
//  OneDimensionalLorentzian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __OneDimensionalLorentzian__
#define __OneDimensionalLorentzian__

#include <iostream>
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"

namespace SpinWaveGenie
{

class OneDimensionalLorentzian : public OneDimensionalShapes
{
public:
  void setFWHM(double InFWHM);
  void setTolerance(double InTolerance) override;
  double getMinimumEnergy() override;
  double getMaximumEnergy() override;
  double getFunction(double frequency, double energy) override;
  std::unique_ptr<OneDimensionalShapes> clone() override;

private:
  double getExponentialFactor();
  double FWHM, Tolerance;
};
}
#endif /* defined(__OneDimensionalLorentzian__) */
