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
#include "Containers/Energies.h"

namespace SpinWaveGenie
{

class EnergyResolutionFunction : public SpinWavePlot{
public:
    EnergyResolutionFunction(){};
    EnergyResolutionFunction(const EnergyResolutionFunction& other);
    EnergyResolutionFunction& operator=(EnergyResolutionFunction& other);
    EnergyResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, Energies energies);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    void setSpinWave (SpinWave SWIn);
    void setResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn);
    const Cell& getCell() const;
    void setEnergies(Energies energies);
    const Energies& getEnergies();
    std::unique_ptr<SpinWavePlot> clone();
    ~EnergyResolutionFunction(){};
private:
    void calculateEnergies();
    std::size_t getBin(double Energy);
    Energies energies;
    std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
    SpinWave SW;
};
}
#endif /* defined(__EnergyResolutionFunction__) */
