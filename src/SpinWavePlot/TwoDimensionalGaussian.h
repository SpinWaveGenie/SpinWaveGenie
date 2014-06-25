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
#include "Genie/SpinWave.h"
#include "SpinWavePlot.h"
#include "Containers/Energies.h"
#include "EnergyResolutionFunction.h"

namespace SpinWaveGenie
{

struct TwoDimGaussian
{
    double a,b,c,tol;
    Vector3 direction;
};

class TwoDimensionResolutionFunction : public SpinWavePlot{
public:
    TwoDimensionResolutionFunction(){};
    TwoDimensionResolutionFunction(TwoDimGaussian& info,SpinWave SW, Energies energies);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    void setTolerance(double toleranceIn);
    std::unique_ptr<SpinWavePlot> clone();
    const Cell& getCell() const;
    const Energies& getEnergies();
    void setEnergies(Energies energies);
    ~TwoDimensionResolutionFunction(){};
private:
    Energies energies;
    double tol;
    double a,b,c;
    double kx,ky,kz;
    Vector3 direction;
    SpinWave SW;
    EnergyResolutionFunction res;
};

}

#endif /* defined(__spin_wave_genie__TwoDimensionalGaussian__) */
