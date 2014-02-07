//
//  OneDimensionalShapes.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __spin_wave_genie__OneDimensionalShapes__
#define __spin_wave_genie__OneDimensionalShapes__

#include <iostream>
#include <memory>

/* Abstract base class */
class OneDimensionalShapes
{
public:
    virtual void setFWHM(double InFWHM) = 0;
    virtual void setTolerance(double InTolerance) = 0;
    virtual double getMinimumEnergy() = 0;
    virtual double getMaximumEnergy() = 0;
    virtual double getFunction(double frequency, double energy) = 0;
    virtual std::unique_ptr<OneDimensionalShapes> clone() = 0;
    virtual ~OneDimensionalShapes(){};
};

#endif /* defined(__spin_wave_genie__OneDimensionalShapes__) */
