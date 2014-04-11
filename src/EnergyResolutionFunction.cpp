//
//  EnergyResolutionFunction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#include "EnergyResolutionFunction.h"

using namespace std;

EnergyResolutionFunction::EnergyResolutionFunction(unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, SpinWave SWIn, double min, double max, size_t points)
{
    //std::cout << "Creating Energy Resolution Function" << std::endl;
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
    //cout << "Energy Points " << EnergyPoints << endl;
    ResolutionFunction = move(ResolutionFunctionIn);
    SW = SWIn;
}

//EnergyResolutionFunction& operator=( const EnergyResolutionFunction& a )
//{
//    ResolutionFunction = move(other.ResolutionFunction0->clone();
//    return *this,
//}

//A& operator=( A&& a )
//{
//    up_ = std::move( a.up_ );
//    return *this,
//}

EnergyResolutionFunction::EnergyResolutionFunction(const EnergyResolutionFunction& other)
{
    std::cout << "Copying Energy Resolution Function" << std::endl;
    MinimumEnergy = other.MinimumEnergy;
    MaximumEnergy = other.MaximumEnergy;
    EnergyPoints = other.EnergyPoints;
    cout << "Energy Points??? " << other.EnergyPoints << endl;
    cout << "Energy Points??? " << EnergyPoints << endl;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
}

EnergyResolutionFunction& EnergyResolutionFunction::operator=(EnergyResolutionFunction other)
{
    std::cout << "Copying Energy Resolution Function" << std::endl;
    MinimumEnergy = other.MinimumEnergy;
    MaximumEnergy = other.MaximumEnergy;
    EnergyPoints = other.EnergyPoints;
    SW = other.SW;
    ResolutionFunction = move(other.ResolutionFunction->clone());
    return *this;
}

std::vector<double> EnergyResolutionFunction::getCut(double kx, double ky, double kz)
{
    //cout << "Energy Points: " << EnergyPoints << endl;
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
        if (isnan(points[k].frequency) || isnan(points[k].intensity))
        {
            
            cout << "found NaN: " << points[k].frequency << " " << points[k].intensity << endl;
        }
        else
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
    }
    
    return fval;
}

double EnergyResolutionFunction::getMinimumEnergy() const
{
    return MinimumEnergy;
}

double EnergyResolutionFunction::getMaximumEnergy() const
{
    return MaximumEnergy;
}

std::size_t EnergyResolutionFunction::getNumberPoints() const
{
    return EnergyPoints;
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


