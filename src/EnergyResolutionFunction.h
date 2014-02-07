//
//  EnergyResolutionFunction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __EnergyResolutionFunction__
#define __EnergyResolutionFunction__

#include <iostream>
#include <memory>
#include "SpinWavePlot.h"
#include "OneDimensionalShapes.h"

class EnergyResolutionFunction : SpinWavePlot{
public:
    EnergyResolutionFunction(){};
    EnergyResolutionFunction(const EnergyResolutionFunction& other) : ResolutionFunction( other.ResolutionFunction->clone() ) {};
    EnergyResolutionFunction& operator=(EnergyResolutionFunction other);
    EnergyResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, double points);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    ~EnergyResolutionFunction(){};
private:
    std::size_t getBin(double Energy);
    double MinimumEnergy,MaximumEnergy;
    std::size_t EnergyPoints;
    SpinWave SW;
    std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
};

#endif /* defined(__EnergyResolutionFunction__) */
