#ifndef __PointsAlongLine__
#define __PointsAlongLine__

#include <iostream>
#include "Positions.h"

class PointsAlongLine
{
public:
    void setFirstPoint(double kx, double ky, double kz);
    void setFinalPoint(double kx, double ky, double kz);
    void setNumberPoints(long points);
    Positions getPoints();
private:
    void calculatePoints();
    Positions Kpoints;
    double kxi,kyi,kzi,kxf,kyf,kzf;
    long numberPoints;    
};

#endif /* defined(__PointsAlongLine__) */
