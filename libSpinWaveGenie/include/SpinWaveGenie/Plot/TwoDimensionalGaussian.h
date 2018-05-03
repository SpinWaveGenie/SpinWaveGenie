//
//  TwoDimensionalGaussian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/21/14.
//
//

#ifndef __spin_wave_genie__TwoDimensionalGaussian__
#define __spin_wave_genie__TwoDimensionalGaussian__

#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include <deque>
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

struct TwoDimGaussian
{
  double a, b, c, tol;
  Eigen::Vector3d direction;
};

class SPINWAVEGENIE_EXPORT TwoDimensionResolutionFunction : public SpinWavePlot
{
public:
  TwoDimensionResolutionFunction() = default;
  TwoDimensionResolutionFunction(const TwoDimGaussian &info, const SpinWave &SW, const Energies &energiesIn);
  TwoDimensionResolutionFunction(const TwoDimensionResolutionFunction & /*other*/) = default;
  std::vector<double> getCut(double kxIn, double kyIn, double kzIn) override;
  void setTolerance(double toleranceIn, int maxEvals = 100000);
  std::unique_ptr<SpinWavePlot> clone() const override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(const Energies &energiesIn) override;

private:
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  Energies energies;
  int maximumEvaluations;
  double tolerance;
  double a, b, c;
  double kx{0.0}, ky{0.0}, kz{0.0};
  Eigen::Vector3d direction;
  EnergyResolutionFunction res;
};
}

#endif /* defined(__spin_wave_genie__TwoDimensionalGaussian__) */
