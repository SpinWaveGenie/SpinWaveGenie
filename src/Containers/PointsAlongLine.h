#ifndef __PointsAlongLine__
#define __PointsAlongLine__

#include <iostream>
#include "Containers/ThreeVectors.h"

class PointsAlongLine
{
public:
    void setFirstPoint(double kx, double ky, double kz);
    void setFinalPoint(double kx, double ky, double kz);
    void setNumberPoints(long points);
    ThreeVectors<double> getPoints();
private:
    void calculatePoints();
    ThreeVectors<double> Kpoints;
    double kxi,kyi,kzi,kxf,kyf,kzf;
    long numberPoints;    
};

#endif /* defined(__PointsAlongLine__) */
