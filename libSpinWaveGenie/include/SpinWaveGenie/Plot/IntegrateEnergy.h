//
//  IntegrateEnergy.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#ifndef __spin_wave_genie__IntegrateEnergy__
#define __spin_wave_genie__IntegrateEnergy__

#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include <deque>
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT IntegrateEnergy : public SpinWavePlot
{
public:
  IntegrateEnergy(const IntegrateEnergy &other);
  IntegrateEnergy(std::unique_ptr<SpinWavePlot> resFunction, const Energies &energies, double delta, double tol = 0.01,
                  int maxEvals = 100000);
  std::vector<double> getCut(double kxIn, double kyIn, double kzIn) override;
  std::unique_ptr<SpinWavePlot> clone() const override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(const Energies &energiesIn) override;

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
