//
//  IntegrateEnergy.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#ifndef __spin_wave_genie__IntegrateEnergy__
#define __spin_wave_genie__IntegrateEnergy__

#include <iostream>
#include <memory>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{
    
    class IntegrateEnergy : public SpinWavePlot
    {
    public:
        IntegrateEnergy(const IntegrateEnergy& other);
        IntegrateEnergy(std::unique_ptr<SpinWavePlot> resFunction, Energies energies, double delta, double tol = 0.01, int maxEval = 100000);
        std::vector<double> getCut(double kx, double ky, double kz);
        //int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
        //static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
        int calculateIntegrand(const int* dim, const double *x,const int* fdim, double *retval);
        static int calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
        std::unique_ptr<SpinWavePlot> clone();
        const Cell& getCell() const;
        const Energies& getEnergies();
        void setEnergies(Energies energies);
        ~IntegrateEnergy(){};
    private:
        std::unique_ptr<SpinWavePlot> resolutionFunction;
        int maximumEvaluations;
        double tolerance,delta;
        double kx,ky,kz;
        Energies centeredEnergies;
    };
}

#endif /* defined(__spin_wave_genie__IntegrateEnergy__) */
