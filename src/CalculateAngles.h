//
//  CalculateAngles.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 4/2/14.
//
//

#ifndef __spin_wave_genie__CalculateAngles__
#define __spin_wave_genie__CalculateAngles__

#include <iostream>
#include <vector>

class CalculateAngles
{
public:
    CalculateAngles(double SA, double SB, double JAB, double JBB, double JBBP, double DB);
    double calculateEnergy(const std::vector<double> &x);
    static double calc(const std::vector<double> &x, std::vector<double> &grad, void * my_func_data);
private:
    double SA,SB,JAB,JBB,JBBP,DB;
    double exchjab(double t0);
    double exchjbb(double t0,double p0,double t1,double p1);
    double exchjbbp(double t0,double p0,double t1,double p1);
};


#endif /* defined(__spin_wave_genie__CalculateAngles__) */
