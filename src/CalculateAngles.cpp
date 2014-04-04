//
//  CalculateAngles.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 4/2/14.
//
//

#include "CalculateAngles.h"
#include <cmath>

CalculateAngles::CalculateAngles(double SAin, double SBin, double JABin, double JBBin, double JBBPin, double DBin)
{
    SA = SAin;
    SB = SBin;
    JAB = JABin;
    JBB = JBBin;
    JBBP = JBBPin;
    DB = DBin;
}

double CalculateAngles::calc(const std::vector<double> &x, std::vector<double> &grad, void * my_func_data)
{
    // Call non-static member function.
    return static_cast<CalculateAngles*>(my_func_data)->calculateEnergy(x);
}

double CalculateAngles::exchjab(double t1)
{
    return -6.0*JAB*SA*SB*cos(t1);
}

double CalculateAngles::exchjbb(double t0,double p0,double t1,double p1)
{
    return -2.0*JBB*SB*SB*(sin(t0)*cos(p0)*sin(t1)*cos(p1)+sin(t0)*sin(p0)*sin(t1)*sin(p1)+cos(t0)*cos(t1));
}

double CalculateAngles::exchjbbp(double t0,double p0,double t1,double p1)
{
    return -2.0*JBBP*SB*SB*(sin(t0)*cos(p0)*sin(t1)*cos(p1)+sin(t0)*sin(p0)*sin(t1)*sin(p1)+cos(t0)*cos(t1));
}

double CalculateAngles::calculateEnergy(const std::vector<double> &angles)
{
    double theta0 = angles[0], theta1 = angles[2], theta2=angles[4], theta3 = angles[6];
    double phi0 = angles[1], phi1 = angles[3], phi2 = angles[5], phi3 = angles[7];
    
    double V0 = 1.0/3.0*DB*SB*SB*pow(-sin(theta0)*cos(phi0)+sin(theta0)*sin(phi0)-cos(theta0),2);
    double V1 = 1.0/3.0*DB*SB*SB*pow(sin(theta1)*cos(phi1)-sin(theta1)*sin(phi1)-cos(theta1),2);
    double V2 = 1.0/3.0*DB*SB*SB*pow(sin(theta2)*cos(phi2)+sin(theta2)*sin(phi2)-cos(theta2),2);
    double V3 = 1.0/3.0*DB*SB*SB*pow(-sin(theta3)*cos(phi3)-sin(theta3)*sin(phi3)-cos(theta3),2);
    
    double JBB02 = exchjbbp(angles[0],angles[1],angles[4],angles[5]);
    double JBB03 = exchjbbp(angles[0],angles[1],angles[6],angles[7]);
    double JBB12 = exchjbbp(angles[2],angles[3],angles[4],angles[5]);
    double JBB13 = exchjbbp(angles[2],angles[3],angles[6],angles[7]);
    
    double JBB01 = exchjbb(angles[0],angles[1],angles[2],angles[3]);
    double JBB23 = exchjbb(angles[4],angles[5],angles[6],angles[7]);
    
    double JAB0 = exchjab(angles[0]);
    double JAB1 = exchjab(angles[2]);
    double JAB2 = exchjab(angles[4]);
    double JAB3 = exchjab(angles[6]);
    
    double energy = V0+V1+V2+V3 + JAB0+JAB1+JAB2+JAB3 + JBB01+JBB02+JBB03+JBB12+JBB13+JBB23;
    //std::cout << "energy = " << energy << std::endl;
    return energy;
}
