//
//  TwoDimensionCut.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//
#include <fstream>
#include "TwoDimensionCut.h"
#include <Eigen/Dense>
#include "EnergyResolutionFunction.h"
#include "Containers/Energies.h"

using std::string; using std::vector; using std::unique_ptr;
using std::cout; using std::endl;

void TwoDimensionCut::setFilename(string name)
{
    Filename = name;
}

void TwoDimensionCut::setPlotObject(unique_ptr<SpinWavePlot> object)
{
    InstrumentResolution = move(object);
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
    figure.setZero(EnergyPoints,Kpoints.size());
    for(auto it = Kpoints.begin(); it != Kpoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();
        
        vector<double> val = InstrumentResolution->getCut(x,y,z);
        size_t m = std::distance(Kpoints.begin(),it);
        cout << m << endl;
        for(int n=0;n<EnergyPoints;n++)
        {
            figure(n,m) = val[n];
        }
    }
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