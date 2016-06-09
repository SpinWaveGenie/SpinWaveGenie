//
//  IntegrateThetaPhi.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#ifndef __spin_wave_genie__IntegrateThetaPhi__
#define __spin_wave_genie__IntegrateThetaPhi__

#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include <deque>
#include <vector>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT IntegrateThetaPhi : public SpinWavePlot
{
public:
  IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tol = 1.0e-4, int maxEvals = 1000)
      : maximumEvaluations(maxEvals), r(0.0), tolerance(tol), resolutionFunction(std::move(object)){};
  IntegrateThetaPhi(const IntegrateThetaPhi &other)
      : maximumEvaluations(other.maximumEvaluations), r(0.0), tolerance(other.tolerance),
        resolutionFunction(other.resolutionFunction->clone()){};
  std::unique_ptr<SpinWavePlot> clone() override;
  const Cell &getCell() const override;
  const Energies &getEnergies() override;
  void setEnergies(const Energies &energies) override;
  std::vector<double> getCut(double kx, double ky, double kz) override;

private:
  std::vector<double> calculateIntegrand(std::deque<double> &x);
  int maximumEvaluations;
  double r, tolerance;
  std::unique_ptr<SpinWavePlot> resolutionFunction;
};
}
#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
