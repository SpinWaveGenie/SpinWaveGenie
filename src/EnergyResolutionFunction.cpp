//
//  EnergyResolutionFunction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#include "EnergyResolutionFunction.h"

using namespace std;

EnergyResolutionFunction::EnergyResolutionFunction(unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, double points)
{
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
    ResolutionFunction = ResolutionFunctionIn->clone();
    SW = SWIn;
}

std::vector<double> EnergyResolutionFunction::getCut(double kx, double ky, double kz)
{
    
    vector<double> fval(EnergyPoints);
    for(int i=0;i!=EnergyPoints;i++)
    {
        fval[i] = 0.0;
    }
    
    SW.createMatrix(kx,ky,kz);
    SW.Calc();
    vector<point> points = SW.getPoints();
    
    for(size_t k=0;k!=points.size();k++)
    {
        double min = points[k].frequency + ResolutionFunction->getMinimumEnergy();
        double max = points[k].frequency + ResolutionFunction->getMaximumEnergy();
        size_t min_bin = getBin(min);
        size_t max_bin = getBin(max);
        
        //cout << min_bin << " " << max_bin << endl;
        for(size_t bin = min_bin; bin <= max_bin; bin++)
        {
            double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)bin/(double)(EnergyPoints-1);
            fval[bin] += points[k].intensity*ResolutionFunction->getFunction(points[k].frequency,energy);
        }
    }
    
    return fval;
}

std::size_t EnergyResolutionFunction::getBin(double Energy)
{
    long bin = ( Energy - MinimumEnergy)*(EnergyPoints-1.0)/(MaximumEnergy-MinimumEnergy);
    
    if (bin < 0)
        bin = 0;
    else if (bin > EnergyPoints)
        bin = EnergyPoints;
    
    return (size_t) bin;
}