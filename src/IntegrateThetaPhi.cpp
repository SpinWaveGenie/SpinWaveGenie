//
//  IntegrateThetaPhi.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#include "IntegrateThetaPhi.h"
#include <vector>
#include "extern/cubature.h"


using namespace std;

IntegrateThetaPhi::IntegrateThetaPhi(EnergyResolutionFunction resFunction, double min, double max, double points, Matrix3 basisVectorsIn)
{
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
    resolutionFunction = resFunction;
    basisVectors = basisVectorsIn;
    tol = 0.001;
}

int IntegrateThetaPhi::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    Vector3 tmp,k;
    double kx,ky,kz;
    double theta = x[0];
    double phi = x[1];
    
    tmp[0] = r*sin(theta)*cos(phi);
    tmp[1] = r*sin(theta)*sin(phi);
    tmp[2] = r*cos(theta);
    
    k = tmp.transpose()*basisVectors.inverse();
    
    //cout << k.transpose() << endl;
    
    vector<double> val = resolutionFunction.getCut(k[0],k[1],k[2]);
    
    for(int i=0;i!=EnergyPoints;i++)
    {
        retval[i] = val[i]/volume;
        //double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
        //cout << energy << " " << val[i] << " ";
    }
    //cout << endl;
    return 0;
}

int IntegrateThetaPhi::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
{
    // Call non-static member function.
    return static_cast<IntegrateThetaPhi*>(data)->calculateIntegrand(dim,x,fdim,retval);
}

std::vector<double> IntegrateThetaPhi::getCut(double kx,double ky, double kz)
{
    std::vector<double> xmin = {0.0,0.0};
    std::vector<double> xmax = {M_PI,2.0*M_PI};
    
    int dim = 2;
    Vector3 dispRLU(kx,ky,kz);
    Vector3 dispAng = dispRLU.transpose()*basisVectors;
    r = dispAng.transpose().norm();
    volume = 4.0*M_PI*r*r;
    
    vector<double> fval(EnergyPoints);
    vector<double> err(EnergyPoints);
    
    hcubature(EnergyPoints,IntegrateThetaPhi::calc, this, dim, &xmin[0], &xmax[0], 0, tol, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    /*for(int i=0;i!=EnergyPoints;i++)
     {
     double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
     cout << energy << " " << fval[i] << " ";
     }
     cout << endl;
     */
    return fval;
    
}
