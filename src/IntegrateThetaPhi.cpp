//
//  IntegrateThetaPhi.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#include "IntegrateThetaPhi.h"
#include <vector>
#include "External/cubature.h"


using namespace std;

IntegrateThetaPhi::IntegrateThetaPhi(double min, double max, double points, Matrix3 basisVectorsIn)
{
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
    basisVectors = basisVectorsIn;
    tol = 0.0001;
}

void IntegrateThetaPhi::setConvolutionObject(EnergyResolutionFunction object)
{
    InstrumentResolution = object;
}


int IntegrateThetaPhi::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    Vector3 tmp,k;
    double theta = x[0];
    double phi = x[1];
    
    //cout << "r= " << r << endl;
    //cout << "theta= " << theta << endl;
    //cout << "phi= " << phi << endl;
    
    tmp[0] = r*sin(theta)*cos(phi);
    tmp[1] = r*sin(theta)*sin(phi);
    tmp[2] = r*cos(theta);
    
    k = basisVectors/(2.0*M_PI)*tmp;
    //cout << tmp.transose() << endl;
    //cout << k.transpose() << endl;
    
    //Vector3 woof = k.transpose()*basisVectors.inverse()*2.0*M_PI;
    //cout << woof.transpose() << endl << endl;
    vector<double> val = InstrumentResolution.getCut(k[0],k[1],k[2]);
    
    //vector<double> values = InstrumentResolution.getCut(0.5,0.0,0.0);
    
    //for (auto it = values.begin(); it!=values.end(); it++)
    //{
    //    cout << (*it) << endl;
    //}
    
    //cout << MinimumEnergy << " " << MaximumEnergy << " " << EnergyPoints << endl;
    double factor = sin(theta)/(4.0*M_PI);
    for(int i=0;i!=EnergyPoints;i++)
    {
        retval[i] = factor*val[i];
        //retval[i] = factor*pow(0.25*sqrt(5.0/M_PI)*(3*cos(theta)*cos(theta)-1.0),1)*0.5*sqrt(3.0/M_PI)*cos(theta)*4.0*M_PI;
        //double S = 1.0;
        //double J = 1.0;
        //double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
        //double frequency = 2.0*J*S*(1.0-cos(tmp[0]));
        //double F = 0.5;
        //cout << retval[i] << " ";
        //cout << factor*(1.0+cos(theta)*cos(theta))*exp(-1.0*(4.0*log(2.0)*pow(energy-frequency,2))/(F*F)) << endl;
        //cout << energy << " " << val[i] << " ";
        //cout << frequency << endl;
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
    Vector3 dispAng = dispRLU.transpose()*2.0*M_PI*basisVectors.inverse();
    r = dispAng.transpose().norm();
    
    //cout << "dispRLU= " << dispRLU.transpose() << endl;
    //cout << "basisVectors= " <<basisVectors << endl;
    //cout << "dispAng = " << dispAng.transpose() << endl;
    
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
