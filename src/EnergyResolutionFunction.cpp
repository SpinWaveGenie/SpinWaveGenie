//
//  EnergyResolutionFunction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//
#include "EnergyResolutionFunction.h"
#include <algorithm>

using namespace std;

EnergyResolutionFunction::EnergyResolutionFunction(unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, size_t points)
{
    //std::cout << "Creating Energy Resolution Function" << std::endl;
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
    this->calculateEnergies();
    //cout << "Energy Points " << EnergyPoints << endl;
    ResolutionFunction = move(ResolutionFunctionIn);
    SW = SWIn;
}

EnergyResolutionFunction::EnergyResolutionFunction(const EnergyResolutionFunction& other)
{
    //std::cout << "Copying Energy Resolution Function" << std::endl;
    MinimumEnergy = other.MinimumEnergy;
    MaximumEnergy = other.MaximumEnergy;
    EnergyPoints = other.EnergyPoints;
    this->calculateEnergies();
    //cout << "Energy Points??? " << other.EnergyPoints << endl;
    //cout << "Energy Points??? " << EnergyPoints << endl;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
}

EnergyResolutionFunction& EnergyResolutionFunction::operator=(EnergyResolutionFunction& other)
{
    std::cout << "Copying Energy Resolution Function" << std::endl;
    MinimumEnergy = other.MinimumEnergy;
    MaximumEnergy = other.MaximumEnergy;
    EnergyPoints = other.EnergyPoints;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
    return *this;
}

void EnergyResolutionFunction::calculateEnergies()
{
    energies.clear();
    energies.reserve(EnergyPoints);
    for (auto bin = 0; bin!=EnergyPoints; bin++)
    {
        energies.push_back(MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)bin/(double)(EnergyPoints-1));
    }
}

std::vector<double> EnergyResolutionFunction::getCut(double kx, double ky, double kz)
{
    //cout << "Energy Points: " << EnergyPoints << endl;
    //cout << MinimumEnergy << " " << MaximumEnergy << endl;
    vector<double> fval(EnergyPoints,0.0);
    
    SW.createMatrix(kx,ky,kz);
    SW.calculate();
    vector<point> points = SW.getPoints();
    
    for(auto pt = points.begin();pt!=points.end();pt++)
    {
        //cout << "k= " << k << endl;
        if (std::isnan(pt->frequency) || std::isnan(pt->intensity))
        {
            
            //cout << "found NaN: " << points[k].frequency << " " << points[k].intensity << endl;
        }
        else
        {
            double min = pt->frequency + ResolutionFunction->getMinimumEnergy();
            double max = pt->frequency + ResolutionFunction->getMaximumEnergy();
            for(size_t index = getBin(min);index!=getBin(max);index++)
            {
                fval[index] += pt->intensity*ResolutionFunction->getFunction(pt->frequency,energies[index]);
            }
        }
    }
    return fval;
}

double EnergyResolutionFunction::getMinimumEnergy() const
{
    return MinimumEnergy;
}

void EnergyResolutionFunction::setMinimumEnergy(double energy)
{
    this->MinimumEnergy = energy;
}

const Cell& EnergyResolutionFunction::getCell() const
{
    return SW.getCell();
}

double EnergyResolutionFunction::getMaximumEnergy() const
{
    return MaximumEnergy;
}

void EnergyResolutionFunction::setMaximumEnergy(double energy)
{
    this->MaximumEnergy = energy;
}

std::size_t EnergyResolutionFunction::getNumberPoints() const
{
    return EnergyPoints;
}

void EnergyResolutionFunction::setNumberPoints(size_t points)
{
    EnergyPoints = points;
}

std::size_t EnergyResolutionFunction::getBin(double Energy)
{
    int bin = round(( Energy - MinimumEnergy)*(EnergyPoints-1.0)/(MaximumEnergy-MinimumEnergy));
    bin = std::max(bin,0);
    return std::min((size_t)bin,EnergyPoints);
}

std::unique_ptr<SpinWavePlot> EnergyResolutionFunction::clone()
{
    return unique_ptr<SpinWavePlot>(new EnergyResolutionFunction(*this));
}


