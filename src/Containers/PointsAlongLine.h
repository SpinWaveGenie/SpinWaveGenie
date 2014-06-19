#ifndef __PointsAlongLine__
#define __PointsAlongLine__

#include <iostream>
#include "Containers/ThreeVectors.h"

namespace SpinWaveGenie
{

//! generates k-points used by the SpinWavePlot routines.
/*!
 This class generates k-points evenly spaced along a line from the first point to the final point.
 The result is stored in a ThreeVectors<double> container.
 */

class PointsAlongLine
{
public:
    //! Set starting point of line. (inclusive)
    //! \param kx H component of k-point,in rlu.
    //! \param ky K component of k-point,in rlu.
    //! \param kz L component of k-point,in rlu.
    void setFirstPoint(double kx, double ky, double kz);
    //! Set final point of line. (inclusive)
    //! \param kx H component of k-point,in rlu.
    //! \param ky K component of k-point,in rlu.
    //! \param kz L component of k-point,in rlu.
    void setFinalPoint(double kx, double ky, double kz);
    //! Set number of k-points along line.
    //! \param points number of k-points.
    void setNumberPoints(long points);
    //! Get resulting k-points along line.
    //! \return k-points along line.
    ThreeVectors<double> getPoints();
private:
    void calculatePoints();
    ThreeVectors<double> Kpoints;
    double kxi,kyi,kzi,kxf,kyf,kzf;
    long numberPoints;    
};
    }
#endif /* defined(__PointsAlongLine__) */
