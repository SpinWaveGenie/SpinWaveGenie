//
//  UniquePositions.cpp
//  positions
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#include "UniquePositions.h"
#include <cmath>

void UniquePositions::insert(double x, double y, double z)
{
    bool unique = true;
    double epsilon = 1.0e-6;
    for(Positions::Iterator it = Positions::begin(); it!= Positions::end(); it++)
    {
        if (std::abs(it->get<0>() - x) < epsilon && std::abs(it->get<1>() - y) < epsilon && std::abs(it->get<2>() - z) < epsilon)
        {
            unique = false;
            break;
        }
    }
    if (unique)
    {
      Positions::insert(x,y,z);
    }
}
