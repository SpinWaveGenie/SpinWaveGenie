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
#include "SpinWaveGenie/Memory.h"
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
  IntegrateAxes(std::unique_ptr<SpinWavePlot> resFunction, HKLDirections directions, double tol = 0.01,
                int maxEval = 100);
  std::vector<double> getCut(double kx, double ky, double kz) override;
  std::unique_ptr<SpinWavePlot> clone() override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(const Energies &energies) override;

private:
  std::unique_ptr<SpinWavePlot> resolutionFunction;
  HKLDirections integrationDirections;
  int maximumEvaluations;
  double tolerance;
  double kx, ky, kz;
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  std::vector<double> xmin, xmax;
};
}
#endif /* defined(__IntegrateAxes__) */
