//
//  IntegrateThetaPhi.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#ifndef __spin_wave_genie__IntegrateThetaPhi__
#define __spin_wave_genie__IntegrateThetaPhi__

#include <iostream>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{

class IntegrateThetaPhi : public SpinWavePlot {
public:
    IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tolerance = 0.01, int maxEvals = 100000);
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
    //int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    //static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    int calculateIntegrand(const int* dim, const double *x,const int* fdim, double *retval);
    static int calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
    int maximumEvaluations;
    double r, tolerance;
    std::unique_ptr<SpinWavePlot> resolutionFunction;
};
}
#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
