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
#include "OneDimensionalShapes.h"

class IntegrateThetaPhi : SpinWavePlot {
    public:
    IntegrateThetaPhi(double min, double max, double points, const Matrix3 basisVectors);
    IntegrateThetaPhi(const IntegrateThetaPhi& other) : InstrumentResolution( other.InstrumentResolution->clone() ) {};
    void setConvolutionObject(std::unique_ptr<OneDimensionalShapes> object);
    std::vector<double> getCut(double kx,double ky, double kz);
    ~IntegrateThetaPhi(){};
    private:
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    double MinimumEnergy,MaximumEnergy;
    double r;
    unsigned EnergyPoints;
    double tol;
    Matrix3 basisVectors;
    std::unique_ptr<OneDimensionalShapes> InstrumentResolution;
};

#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
