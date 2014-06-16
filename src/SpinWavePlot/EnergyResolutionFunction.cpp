//
//  EnergyResolutionFunction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//
#include "EnergyResolutionFunction.h"
#include "Containers/Results.h"
#include <algorithm>

using namespace std;

EnergyResolutionFunction::EnergyResolutionFunction(unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, Energies energiesIn)
{
    //std::cout << "Creating Energy Resolution Function" << std::endl;
    this->energies = energiesIn;
    //cout << "Energy Points " << EnergyPoints << endl;
    ResolutionFunction = move(ResolutionFunctionIn);
    SW = SWIn;
}

EnergyResolutionFunction::EnergyResolutionFunction(const EnergyResolutionFunction& other)
{
    //std::cout << "Copying Energy Resolution Function" << std::endl;
    energies = other.energies;
    //cout << "Energy Points??? " << other.EnergyPoints << endl;
    //cout << "Energy Points??? " << EnergyPoints << endl;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
}

EnergyResolutionFunction& EnergyResolutionFunction::operator=(EnergyResolutionFunction& other)
{
    std::cout << "Copying Energy Resolution Function" << std::endl;
    energies = other.energies;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
    return *this;
}


std::vector<double> EnergyResolutionFunction::getCut(double kx, double ky, double kz)
{
    //cout << "Energy Points: " << EnergyPoints << endl;
    //cout << MinimumEnergy << " " << MaximumEnergy << endl;
    size_t EnergyPoints = energies.size();
    vector<double> fval(EnergyPoints,0.0);
    
    SW.createMatrix(kx,ky,kz);
    SW.calculate();
    Results points = SW.getPoints();
    
    for(auto pt = points.begin();pt!=points.end();pt++)
    {
        if (std::isnan(pt->frequency) || std::isnan(pt->intensity))
        {
            cout << "found NaN: " << pt->frequency << " " << pt->intensity << endl;
        }
        else
        {
            double min = pt->frequency + ResolutionFunction->getMinimumEnergy();
            double max = pt->frequency + ResolutionFunction->getMaximumEnergy();
            size_t UpperBound = energies.getUpperBound(max);
            //cout << min << " " << energies.getLowerBound(min) << " " << max << " " << UpperBound << endl;
            for(size_t index = energies.getLowerBound(min);index!=UpperBound;index++)
            {
                fval[index] += pt->intensity*ResolutionFunction->getFunction(pt->frequency,energies[index]);
            }
        }
    }
    return fval;
}

const Cell& EnergyResolutionFunction::getCell() const
{
    return SW.getCell();
}

const Energies& EnergyResolutionFunction::getEnergies()
{
    return energies;
}

void EnergyResolutionFunction::setEnergies(Energies energiesIn)
{
    energies = energiesIn;
}

std::unique_ptr<SpinWavePlot> EnergyResolutionFunction::clone()
{
    return unique_ptr<SpinWavePlot>(new EnergyResolutionFunction(*this));
}


