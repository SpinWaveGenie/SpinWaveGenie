//
//  HKLDirection.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/22/14.
//
//

#include "SpinWaveGenie/Containers/HKLDirections.h"

namespace SpinWaveGenie
{

void HKLDirections::addDirection(double v0, double v1, double v2, double delta)
{
    Axis direction;
    direction.v0 = v0;
    direction.v1 = v1;
    direction.v2 = v2;
    direction.delta = delta;
    
    integrationDirections.push_back(direction);
}

void HKLDirections::addDirection(int direction, double delta)
{
    double v0 = 0.0, v1=0.0, v2=0.0;
    switch(direction)
    {
        case 0:
            v0 = 1.0;
            break;
        case 1:
            v1 = 1.0;
            break;
        case 2:
            v2 = 1.0;
            break;
    }
    
    this->addDirection(v0, v1, v2, delta);
}

}