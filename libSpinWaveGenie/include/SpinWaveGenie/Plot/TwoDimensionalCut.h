//
//  TwoDimensionalCut.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//

#ifndef __spin_wave_genie__TwoDimensionalCut__
#define __spin_wave_genie__TwoDimensionalCut__

#include <memory>
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"

namespace SpinWaveGenie
{
    
class TwoDimensionalCut
{
public:
    TwoDimensionalCut();
    TwoDimensionalCut(const TwoDimensionalCut& Other);
    void setFilename(std::string name);
    void setPoints(ThreeVectors<double> pos);
    void setEnergyPoints(double min, double max, size_t numberpoints);
    void setPlotObject(std::unique_ptr<SpinWavePlot> object);
    Eigen::MatrixXd getMatrix();
    void save();
    ~TwoDimensionalCut();
private:
    class CutImpl;
    std::unique_ptr<CutImpl> m_p;
};

}
#endif /* defined(__spin_wave_genie__TwoDimensionalCut__) */
