//
//  TwoDimensionCut.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//
#include <fstream>
#include <Eigen/Dense>
#include "SpinWaveGenie/SpinWavePlot/TwoDimensionCut.h"
#include "SpinWaveGenie/SpinWavePlot/EnergyResolutionFunction.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "External/ezRateProgressBar.hpp"

namespace SpinWaveGenie
{

using std::string; using std::vector; using std::unique_ptr;
using std::cout; using std::endl;

void TwoDimensionCut::setFilename(string name)
{
    Filename = name;
}

void TwoDimensionCut::setPlotObject(unique_ptr<SpinWavePlot> object)
{
    InstrumentResolution = move(object);
    EnergyPoints = InstrumentResolution->getEnergies().size();
}

void TwoDimensionCut::setPoints(ThreeVectors<double> pts)
{
    Kpoints.clear();
    for(auto it = pts.begin(); it!= pts.end(); it++)
    {
        Kpoints.insert(it->get<0>(),it->get<1>(),it->get<2>());
    }
}

void TwoDimensionCut::setEnergyPoints(double min, double max, size_t points)
{
    InstrumentResolution->setEnergies(Energies(min, max, points));
    EnergyPoints = points;
}

Eigen::MatrixXd TwoDimensionCut::getMatrix()
{
    Eigen::MatrixXd figure;
    figure.setZero(Kpoints.size(),EnergyPoints);
    ez::ezRateProgressBar<int> p(Kpoints.size());
    p.units = "Q-points";
    p.start();
    for(auto it = Kpoints.begin(); it != Kpoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();
        
        vector<double> val = InstrumentResolution->getCut(x,y,z);
        size_t m = std::distance(Kpoints.begin(),it);
        p.update(m);
        for(int n=0;n<EnergyPoints;n++)
        {
            figure(m,n) = val[n];
        }
    }
    p.update(Kpoints.size());
    return figure;
}

void TwoDimensionCut::save()
{
    Eigen::MatrixXd figure = this->getMatrix();
    std::ofstream file(Filename);
    if (file.is_open())
    {
        file << figure;
    }
    file << endl;
    file.close();
}

}