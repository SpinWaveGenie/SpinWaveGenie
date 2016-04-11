//
//  IntegrateEnergy.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#ifndef __spin_wave_genie__IntegrateEnergy__
#define __spin_wave_genie__IntegrateEnergy__

#include <iostream>
#include "SpinWaveGenie/Memory.h"
#include <vector>
#include <deque>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{

class IntegrateEnergy : public SpinWavePlot
{
public:
  IntegrateEnergy(const IntegrateEnergy &other);
  IntegrateEnergy(std::unique_ptr<SpinWavePlot> resFunction, const Energies &energies, double delta, double tol = 0.01,
                  int maxEval = 100000);
  std::vector<double> getCut(double kx, double ky, double kz) override;
  std::unique_ptr<SpinWavePlot> clone() override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(Energies energies) override;

private:
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  std::unique_ptr<SpinWavePlot> resolutionFunction;
  int maximumEvaluations;
  double tolerance, delta;
  double kx, ky, kz;
  Energies centeredEnergies;
};
}

#endif /* defined(__spin_wave_genie__IntegrateEnergy__) */
