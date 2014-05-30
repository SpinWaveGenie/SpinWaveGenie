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
#include "SpinWavePlot.h"
#include "Containers/HKLDirections.h"
#include "Containers/Energies.h"

class IntegrateAxes : public SpinWavePlot
{
public:
    IntegrateAxes(const IntegrateAxes& other);
    IntegrateAxes(std::unique_ptr<SpinWavePlot> resFunction, HKLDirections directions, double tol);
    std::vector<double> getCut(double kx, double ky, double kz);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::unique_ptr<SpinWavePlot> clone();
    const Cell& getCell() const;
    Energies getEnergies();
    void setEnergies(Energies energies);
    ~IntegrateAxes(){};
private:
    std::unique_ptr<SpinWavePlot> resolutionFunction;
    HKLDirections integrationDirections;
    double tolerance, volume;
    double kx,ky,kz;
};

#endif /* defined(__IntegrateAxes__) */
