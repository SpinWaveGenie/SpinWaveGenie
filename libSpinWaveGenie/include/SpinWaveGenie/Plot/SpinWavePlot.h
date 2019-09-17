#ifndef __SpinWavePlot__
#define __SpinWavePlot__

#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Plot/OneDimensionalGaussian.h"
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

/* Abstract base class */
class SPINWAVEGENIE_EXPORT SpinWavePlot
{
public:
  virtual std::unique_ptr<SpinWavePlot> clone() const = 0;
  virtual const Cell &getCell() const = 0;
  virtual const Energies &getEnergies() = 0;
  virtual void setEnergies(const Energies &energies) = 0;
  virtual std::vector<double> getCut(double kx, double ky, double kz) = 0; // returns constant-Q cut
  virtual ~SpinWavePlot() = default;
};
} // namespace SpinWaveGenie
#endif /* defined(__SpinWavePlot__) */
