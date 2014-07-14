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
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/SpinWavePlot/SpinWavePlot.h"
#include "SpinWaveGenie/SpinWavePlot/OneDimensionalShapes.h"

namespace SpinWaveGenie
{

class TwoDimensionCut
{
public:
    void setFilename(std::string name);
    void setPoints(ThreeVectors<double> pos);
    void setEnergyPoints(double min, double max, size_t numberpoints);
    void setPlotObject(std::unique_ptr<SpinWavePlot> object);
    Eigen::MatrixXd getMatrix();
    void save();
private:
    size_t EnergyPoints;
    std::unique_ptr<SpinWavePlot> InstrumentResolution;
    std::string Filename;
    ThreeVectors<double> Kpoints;
};

}
#endif /* defined(__spin_wave_genie__TwoDimensionCut__) */
