//
//  IntegrateThetaPhi.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#ifndef __spin_wave_genie__IntegrateThetaPhi__
#define __spin_wave_genie__IntegrateThetaPhi__

#include <deque>
#include <vector>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{

class IntegrateThetaPhi : public SpinWavePlot {
public:
    IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tolerance = 1.0e-4, int maxEvals = 1000);
    IntegrateThetaPhi(const IntegrateThetaPhi& other) : maximumEvaluations(other.maximumEvaluations),
                                                        tolerance(other.tolerance),
                                                        resolutionFunction(move(other.resolutionFunction->clone())) {};
    std::unique_ptr<SpinWavePlot> clone();
    const Cell& getCell() const;
    const Energies& getEnergies();
    void setEnergies(Energies energies);
    std::vector<double> getCut(double kx,double ky, double kz);
    ~IntegrateThetaPhi(){};
private:
    std::vector<double> calculateIntegrand(std::deque<double>& x);
    int maximumEvaluations;
    double r, tolerance;
    std::unique_ptr<SpinWavePlot> resolutionFunction;
};
}
#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
