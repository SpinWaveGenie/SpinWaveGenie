//
//  File.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/30/14.
//
//

#ifndef __OneDimensionalGaussian__
#define __OneDimensionalGaussian__

#include <iostream>

class OneDimensionalGaussian
{
public:
    void setFWHM(double InFWHM);
    void setTolerance(double InTolerance);
    double getMinimumEnergy();
    double getMaximumEnergy();
    double getFunction(double frequency, double energy);
private:
    double getExponentialFactor();
    double FWHM,Tolerance;
};

#endif /* defined(__OneDimensionalGaussian__) */
