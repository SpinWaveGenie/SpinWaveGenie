//
//  Positions.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#include "Positions.h"
#include <cmath>

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