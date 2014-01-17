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

using std::string; using std::vector;
using std::cout; using std::endl;

void TwoDimensionCut::setFilename(string name)
{
    Filename = name;
}

void TwoDimensionCut::setConvolutionObject(OneDimGaussian object)
{
    InstrumentResolution = object;
}

void TwoDimensionCut::setPoints(Positions pts)
{
    for( Positions::Iterator it = pts.begin(); it!= pts.end(); it++)
    {
        Kpoints.insert(it->get<0>(),it->get<1>(),it->get<2>());
    }
}

void TwoDimensionCut::setEnergyPoints(double min, double max, size_t points)
{
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
}

void TwoDimensionCut::save()
{
    Eigen::MatrixXd figure;
    figure.setZero(EnergyPoints,Kpoints.size());
    EnergyResolutionFunction scan(InstrumentResolution,MinimumEnergy,MaximumEnergy,EnergyPoints);
    
    for(Positions::Iterator it = Kpoints.begin(); it != Kpoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();

        vector<double> val = scan.getCut(x,y,z);
        size_t m = std::distance(Kpoints.begin(),it);
        for(int n=0;n<EnergyPoints;n++)
        {
            figure(n,m) = val[n];
        }
    }
    std::ofstream file(Filename);
    if (file.is_open())
    {
        file << figure;
    }
    file << endl;
}