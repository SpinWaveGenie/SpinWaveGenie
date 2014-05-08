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

struct axes_info
{
    bool x,y,z;
    double dx,dy,dz;
    Vector3 v1,v2,v3;
    double s1,s2,s3;
    double tol;
};

class IntegrateAxes
{
public:
    IntegrateAxes(axes_info info, std::unique_ptr<SpinWavePlot> resFunction);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::vector<double> getCut(double kx, double ky, double kz);
    ~IntegrateAxes(){};
private:
    axes_info info;
    double kx,ky,kz;
    double volume;
    std::size_t EnergyPoints;
    double MinimumEnergy,MaximumEnergy;
    std::unique_ptr<SpinWavePlot> ResolutionFunction;
};

#endif /* defined(__IntegrateAxes__) */
