//
//  Energies.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/27/14.
//
//

#ifndef __spin_wave_genie__Energies__
#define __spin_wave_genie__Energies__

#include <iostream>
#include <vector>

class Energies
{
public:
    Energies() {};
    Energies(double minimum, double maximum, std::size_t numberPoints);
    void insert(double energy);
    std::size_t getLowerBound(double energy);
    std::size_t getUpperBound(double energy);
    std::size_t size();
    void clear();
    double* data();
    const double& operator[](std::size_t position);
    typedef std::vector<double>::iterator Iterator;
    Iterator begin();
    Iterator end();
private:
    std::vector<double> energies;
};

inline std::size_t Energies::size()
{
    return energies.size();
}

inline const double& Energies::operator[](std::size_t bin)
{
    return energies[bin];
}

inline Energies::Iterator Energies::begin()
{
    return energies.begin();
}

inline Energies::Iterator Energies::end()
{
    return energies.end();
}

#endif /* defined(__spin_wave_genie__Energies__) */
