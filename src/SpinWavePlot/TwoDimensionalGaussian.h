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

namespace SpinWaveGenie
{

struct TwoDimGaussian
{
    double a,b,c;
    unsigned direction;
    double tol;
    SpinWave SW;
};

class TwoDimensionResolutionFunction : public SpinWavePlot{
public:
    TwoDimensionResolutionFunction(){};
    TwoDimensionResolutionFunction(TwoDimGaussian& info, double min, double max, double points);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    ~TwoDimensionResolutionFunction(){};
private:
    double MinimumEnergy,MaximumEnergy,tol;
    double a,b,c;
    double kx,ky,kz;
    unsigned EnergyPoints,direction;
    SpinWave SW;
};

}

#endif /* defined(__spin_wave_genie__TwoDimensionalGaussian__) */
