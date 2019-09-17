//
//  IntegrateAxes.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/7/14.
//
//

#ifndef __IntegrateAxes__
#define __IntegrateAxes__

#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Containers/HKLDirections.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include <deque>
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT IntegrateAxes : public SpinWavePlot
{
public:
  IntegrateAxes(const IntegrateAxes &other);
  IntegrateAxes(std::unique_ptr<SpinWavePlot> &&resFunction, const HKLDirections &directions, double tol = 0.01,
                int maxEvals = 100);
  std::vector<double> getCut(double kxIn, double kyIn, double kzIn) override;
  std::unique_ptr<SpinWavePlot> clone() const override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(const Energies &energiesIn) override;

private:
  std::unique_ptr<SpinWavePlot> resolutionFunction;
  HKLDirections integrationDirections;
  int maximumEvaluations;
  double tolerance;
  double kx{0.0}, ky{0.0}, kz{0.0};
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  std::vector<double> xmin, xmax;
};
} // namespace SpinWaveGenie
#endif /* defined(__IntegrateAxes__) */
