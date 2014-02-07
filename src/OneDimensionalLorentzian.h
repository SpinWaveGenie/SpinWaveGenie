//
//  OneDimensionalLorentzian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __OneDimensionalLorentzian__
#define __OneDimensionalLorentzian__

#include <iostream>
#include <memory>
#include "OneDimensionalShapes.h"

class OneDimensionalLorentzian: public OneDimensionalShapes
{
public:
    void setFWHM(double InFWHM);
    void setTolerance(double InTolerance);
    double getMinimumEnergy();
    double getMaximumEnergy();
    double getFunction(double frequency, double energy);
    std::unique_ptr<OneDimensionalShapes> clone();
    ~OneDimensionalLorentzian(){};
private:
    double getExponentialFactor();
    double FWHM,Tolerance;
};

#endif /* defined(__OneDimensionalLorentzian__) */
