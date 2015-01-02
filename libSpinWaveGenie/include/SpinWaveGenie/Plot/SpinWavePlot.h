#ifndef __SpinWavePlot__
#define __SpinWavePlot__

#include <iostream>
#include <string>
#include <vector>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{

/* Abstract base class */
class SpinWavePlot
{
public:
  virtual std::unique_ptr<SpinWavePlot> clone() = 0;
  virtual const Cell &getCell() const = 0;
  virtual const Energies &getEnergies() = 0;
  virtual void setEnergies(Energies energies) = 0;
  virtual std::vector<double> getCut(double kx, double ky, double kz) = 0; // returns constant-Q cut
  virtual ~SpinWavePlot(){};
};
}
#endif /* defined(__SpinWavePlot__) */
