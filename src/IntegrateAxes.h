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

struct Axis
{
    double v0,v1,v2,delta;
};


class IntegrateAxes
{
public:
    IntegrateAxes(std::unique_ptr<SpinWavePlot> resFunction, double tol);
    void addDirection(int direction, double delta);
    void addDirection(double v0, double v1, double v2, double delta);
    std::vector<double> getCut(double kx, double ky, double kz);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    ~IntegrateAxes(){};
private:
    std::unique_ptr<SpinWavePlot> resolutionFunction;
    std::vector<Axis> integrationDirections;
    double tolerance, volume, minimumEnergy, maximumEnergy;
    std::size_t energyPoints;
    double kx,ky,kz;
    
};

#endif /* defined(__IntegrateAxes__) */
