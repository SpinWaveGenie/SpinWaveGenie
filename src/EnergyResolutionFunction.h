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
    EnergyResolutionFunction& operator=(EnergyResolutionFunction& other);
    EnergyResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, std::size_t points);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    double getMinimumEnergy() const;
    void setMinimumEnergy(double energy);
    double getMaximumEnergy() const;
    void setMaximumEnergy(double energy);
    std::size_t getNumberPoints() const;
    void setNumberPoints(std::size_t points);
    const Cell& getCell() const;
    std::unique_ptr<SpinWavePlot> clone();
    ~EnergyResolutionFunction(){};
private:
    void calculateEnergies();
    std::size_t getBin(double Energy);
    std::size_t EnergyPoints;
    std::vector<double> energies;
    double MinimumEnergy,MaximumEnergy;
    std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
    SpinWave SW;
};

#endif /* defined(__EnergyResolutionFunction__) */
