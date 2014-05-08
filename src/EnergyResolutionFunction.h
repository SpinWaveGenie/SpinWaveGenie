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
#include "Cell/Cell.h"
#include "OneDimensionalShapes.h"

class EnergyResolutionFunction : public SpinWavePlot{
public:
    EnergyResolutionFunction(){};
    EnergyResolutionFunction(const EnergyResolutionFunction& other);
    EnergyResolutionFunction& operator=(EnergyResolutionFunction other);
    EnergyResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, std::size_t points);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    double getMinimumEnergy() const;
    double getMaximumEnergy() const;
    const Cell& getCell() const;
    std::size_t getNumberPoints() const;
    ~EnergyResolutionFunction(){};
private:
    std::size_t getBin(double Energy);
    std::size_t EnergyPoints;
    double MinimumEnergy,MaximumEnergy;
    std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
    SpinWave SW;
};

#endif /* defined(__EnergyResolutionFunction__) */
