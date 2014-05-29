//
//  Energies.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/27/14.
//
//

#include "Energies.h"
#include <algorithm>


Energies::Energies(double minimum, double maximum, std::size_t numberPoints)
{
    energies.reserve(numberPoints);
    for (auto bin = 0; bin!=numberPoints; bin++)
    {
        energies.push_back(minimum + (maximum - minimum)*(double)bin/(double)(numberPoints-1));
    }
}

void Energies::insert(double energy)
{
    energies.push_back(energy);
    std::sort(energies.begin(),energies.end());
}

std::size_t Energies::getLowerBound(double energy)
{
    return std::distance(energies.begin(),std::lower_bound(energies.begin(),energies.end(),energy));
}

std::size_t Energies::getUpperBound(double energy)
{
    return std::distance(energies.begin(),std::upper_bound(energies.begin(),energies.end(),energy));
}

void Energies::clear()
{
    energies.clear();
}
