//
//  TwoDimensionalCut.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//

#ifndef __spin_wave_genie__TwoDimensionalCut__
#define __spin_wave_genie__TwoDimensionalCut__

#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Export.h"

#include <memory>

namespace SpinWaveGenie
{

class SpinWavePlot;

class SPINWAVEGENIE_EXPORT TwoDimensionalCut
{
public:
  void setFilename(const std::string &name);
  void setPoints(const ThreeVectors<double> &pts);
  void setEnergyPoints(double min, double max, size_t points);
  void setPlotObject(std::unique_ptr<SpinWavePlot> object);
  Eigen::MatrixXd getMatrix();
  void save();
private:
  std::string filename;
  std::unique_ptr<SpinWaveGenie::SpinWavePlot> cut;
  SpinWaveGenie::ThreeVectors<double> points;
};
}
#endif /* defined(__spin_wave_genie__TwoDimensionalCut__) */
