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
#include "Containers/Positions.h"
#include "Genie/SpinWave.h"
#include "SpinWavePlot.h"

class TwoDimensionCut
{
public:
    void setFilename(std::string name);
    void setPoints(Positions pos);
    void setGenie(SpinWave SW);
    void setConvolutionObject(OneDimGaussian object);
    void save();
private:
    double MaximumEnergy,MinimumEnergy;
    size_t EnergyPoints;
    OneDimGaussian InstrumentResolution;
    SpinWave Genie;
    std::string Filename;
    Positions Kpoints;
};


#endif /* defined(__spin_wave_genie__TwoDimensionCut__) */
