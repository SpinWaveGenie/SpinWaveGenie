//
//  TwoDimensionalGaussian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/21/14.
//
//

#ifndef __spin_wave_genie__TwoDimensionalGaussian__
#define __spin_wave_genie__TwoDimensionalGaussian__

#include <iostream>
#include <vector>
#include <deque>
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"

namespace SpinWaveGenie
{

struct TwoDimGaussian
{
  double a, b, c, tol;
  Vector3 direction;
};

class TwoDimensionResolutionFunction : public SpinWavePlot
{
public:
  TwoDimensionResolutionFunction(){};
  TwoDimensionResolutionFunction(TwoDimGaussian &info, SpinWave SW, Energies energies);
  TwoDimensionResolutionFunction(const TwoDimensionResolutionFunction &other) = default;
  std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
  void setTolerance(double tol, int maxEvals = 100000);
  std::unique_ptr<SpinWavePlot> clone();
  const Cell &getCell() const;
  const Energies &getEnergies();
  void setEnergies(Energies energies);
  ~TwoDimensionResolutionFunction(){};

private:
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  Energies energies;
  int maximumEvaluations;
  double tolerance;
  double a, b, c;
  double kx, ky, kz;
  Vector3 direction;
  EnergyResolutionFunction res;
  double xmax;
};
}

#endif /* defined(__spin_wave_genie__TwoDimensionalGaussian__) */
