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
#include "SpinWavePlot.h"

class IntegrateThetaPhi : SpinWavePlot {
    public:
    IntegrateThetaPhi( EnergyResolutionFunction resFunction, double min, double max, double points, const Matrix3 basisVectors );
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::vector<double> getCut(double kx,double ky, double kz);
    ~IntegrateThetaPhi(){};
    private:
    double MinimumEnergy,MaximumEnergy;
    double r;
    unsigned EnergyPoints;
    double tol,volume;
    EnergyResolutionFunction resolutionFunction;
    Matrix3 basisVectors;
};

#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
