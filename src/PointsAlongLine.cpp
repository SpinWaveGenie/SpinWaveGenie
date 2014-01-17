#include "PointsAlongLine.h"

void PointsAlongLine::setFirstPoint(double kx, double ky, double kz)
{
    kxi = kx;
    kyi = ky;
    kzi = kz;
}

void PointsAlongLine::setFinalPoint(double kx, double ky, double kz)
{
    kxf = kx;
    kyf = ky;
    kzf = kz;
}

void PointsAlongLine::setNumberPoints(long points)
{
    numberPoints = points;
}

void PointsAlongLine::calculatePoints()
{
    for(int m=0;m<numberPoints;m++)
    {
        if(numberPoints==1)
        {
            Kpoints.insert(kxi,kyi,kzi);
            //std::cout << kxi << " " << kyi << " " << kzi << std::endl;
        }
        else
        {
            double x = kxi + (kxf-kxi)*m/(numberPoints-1);
            double y = kyi + (kyf-kyi)*m/(numberPoints-1);
            double z = kzi + (kzf-kzi)*m/(numberPoints-1);
            Kpoints.insert(x,y,z);
            //std::cout << x << " " << y << " " << z << std::endl;
        }
    }
}



Positions PointsAlongLine::getPoints()
{
    calculatePoints();
    return Kpoints;
}