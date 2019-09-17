//
//  OneDimensionalPseudoVoigt.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/6/14.
//
//

#ifndef __OneDimensionalPseudoVoigt__
#define __OneDimensionalPseudoVoigt__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/OneDimensionalLorentzian.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT OneDimensionalPseudoVoigt : public OneDimensionalShapes
{
public:
  OneDimensionalPseudoVoigt();
  OneDimensionalPseudoVoigt(const OneDimensionalPseudoVoigt &other)
      : eta(other.eta), tolerance(other.tolerance), Lorentzian(other.Lorentzian->clone()),
        Gaussian(other.Gaussian->clone()){};
  void setEta(double InEta);
  void setFWHM(double InFWHM);
  void setTolerance(double InTolerance) override;
  double getMinimumEnergy() override;
  double getMaximumEnergy() override;
  void setFrequency(double frequency) override;
  double getFunction(double energy) override;
  std::unique_ptr<OneDimensionalShapes> clone() const override;

private:
  double eta{0.0}, tolerance{0.001};
  std::unique_ptr<OneDimensionalShapes> Lorentzian{std::make_unique<OneDimensionalGaussian>()},
      Gaussian{std::make_unique<OneDimensionalLorentzian>()};
};
} // namespace SpinWaveGenie
#endif /* defined(__OneDimensionalPseudoVoigt__) */
