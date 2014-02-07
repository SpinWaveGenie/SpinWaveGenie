//
//  TwoDimensionCut.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//

#ifndef __spin_wave_genie__TwoDimensionCut__
#define __spin_wave_genie__TwoDimensionCut__

#include <iostream>
#include <memory>
#include "Containers/Positions.h"
#include "Genie/SpinWave.h"
#include "SpinWavePlot.h"
#include "OneDimensionalShapes.h"

class TwoDimensionCut
{
public:
    
    void setFilename(std::string name);
    void setPoints(Positions pos);
    void setEnergyPoints(double min, double max, size_t numberpoints);
    void setConvolutionObject(std::unique_ptr<OneDimensionalShapes> object);
    void setSpinWave(SpinWave SWIn);
    void save();
private:
    double MaximumEnergy,MinimumEnergy;
    size_t EnergyPoints;
    std::unique_ptr<OneDimensionalShapes> InstrumentResolution;
    SpinWave SW;
    std::string Filename;
    Positions Kpoints;
};


#endif /* defined(__spin_wave_genie__TwoDimensionCut__) */
