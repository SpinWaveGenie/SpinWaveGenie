//
//  Positions.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#include "Positions.h"
#include <cmath>

bool Positions::operator==(Positions& other)
{
    for(auto MyPosition = this->begin(); MyPosition!= this->end(); MyPosition++)
    {
        bool ContainsPosition = false;
        for(auto OtherPosition = other.begin(); OtherPosition!= other.end(); OtherPosition++)
        {
            double x = MyPosition->get<0>() - OtherPosition->get<0>();
            double y = MyPosition->get<1>() - OtherPosition->get<1>();
            double z = MyPosition->get<2>() - OtherPosition->get<2>();
            double norm = sqrt(x*x+y*y+z*z);
            if (norm < 1.0e-8)
                ContainsPosition = true;
        }
        if(ContainsPosition == false)
            return ContainsPosition;
    }
    return true;
}


void Positions::insert(double x, double y, double z)
{
  valuesX.push_back(x);
  valuesY.push_back(y);
  valuesZ.push_back(z);
}

size_t Positions::size()
{
    return valuesX.size();
}

Positions::Iterator Positions::begin()
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.begin(),valuesY.begin(),valuesZ.begin()));
}

Positions::Iterator Positions::end()
{
  return boost::make_zip_iterator(boost::make_tuple(valuesX.end(),valuesY.end(),valuesZ.end()));
}