//
//  IntegrateAxes.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/7/14.
//
//

#ifndef __IntegrateAxes__
#define __IntegrateAxes__

#include <iostream>
#include <memory>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/HKLDirections.h"
#include "SpinWaveGenie/Containers/Energies.h"

namespace SpinWaveGenie
{

class IntegrateAxes : public SpinWavePlot
{
public:
    IntegrateAxes(const IntegrateAxes& other);
    IntegrateAxes(std::unique_ptr<SpinWavePlot> resFunction, HKLDirections directions, double tol = 0.01, int maxEval = 100000);
    std::vector<double> getCut(double kx, double ky, double kz);
    int calculateIntegrand(const int* dim, const double *x,const int* fdim, double *retval);
    static int calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);
    std::unique_ptr<SpinWavePlot> clone();
    const Cell& getCell() const;
    const Energies& getEnergies();
    void setEnergies(Energies energies);
    ~IntegrateAxes(){};
private:
    std::unique_ptr<SpinWavePlot> resolutionFunction;
    HKLDirections integrationDirections;
    int maximumEvaluations;
    double tolerance,volume;
    double kx,ky,kz;
    std::vector<double> xmin,xmax;
};
}
#endif /* defined(__IntegrateAxes__) */
