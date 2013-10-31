//
//  SpinWavePlot.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 10/30/13.
//
//
#include <cmath>
#include <functional>
#include "SpinWavePlot.h"
#include "SpinWave.h"
#include "extern/cubature.h"

using namespace std;

TwoDimensionResolutionFunction::TwoDimensionResolutionFunction(SW_Builder& builderInput, double min, double max, double points)
{
    builder = builderInput;
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
}

int TwoDimensionResolutionFunction::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *fval)
{
    assert(fdim == EnergyPoints);
    assert(dim == 1);
    
    double a,b,c;
    a = 1109.0;
    b = 0.0;
    c = 0.48;
    
    for(int i=0;i!=EnergyPoints;i++)
    {
        fval[i] = 0.0;
    }

    SpinWave* test = builder.Create_Element(kx,x[0],kz);
    test->Calc();
    vector<double> frequencies = test->Get_Frequencies();
    vector<double> intensities = test->Get_Intensities();
    double sigma_energy = 1.0/sqrt(2.0*c);

    double u = x[0] - ky;
    for(size_t k=0;k!=frequencies.size();k++)
    {
        //cout << frequencies[k] << " " << intensities[k] << "  ";
        long min_bin = (long) (frequencies[k]-7.0*sigma_energy - MinimumEnergy)*(EnergyPoints-1)/(MaximumEnergy-MinimumEnergy);
        long max_bin = (long) (frequencies[k]+7.0*sigma_energy - MinimumEnergy)*(EnergyPoints-1)/(MaximumEnergy-MinimumEnergy);
        
        if (min_bin < 0)
            min_bin = 0;
        if (max_bin > EnergyPoints)
            max_bin = EnergyPoints;
                
        for(int i=min_bin;i!=max_bin;i++)
        {
            double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
            fval[i] += intensities[k]*exp(-c*pow(frequencies[k]-energy,2))*exp(-2.0*b*(frequencies[k]-energy)*u)*exp(-a*pow(u,2));
        }
    }
    //cout << endl;
    /*for(int i=0;i!=EnergyPoints;i++)
    {
        double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
        cout << energy << " " << fval[i] << " ";
    }
    cout << endl;*/
    
    return 0;
}

int TwoDimensionResolutionFunction::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
{
    // Call non-static member function.
    return static_cast<TwoDimensionResolutionFunction*>(data)->calculateIntegrand(dim,x,fdim,retval);
}

std::vector<double> TwoDimensionResolutionFunction::getCut(double kxIn, double kyIn, double kzIn)
{
    kx = kxIn;
    ky = kyIn;
    kz = kzIn;
    
    double xmin, xmax;
    double tol = 1.0e-4;
    
    xmin = ky - 0.2;
    xmax = ky + 0.2;
    
    vector<double> fval(EnergyPoints);
    vector<double> err(EnergyPoints);
    
    hcubature(EnergyPoints,TwoDimensionResolutionFunction::calc, this, 1, &xmin, &xmax, 0, tol, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    return fval;
}


IntegrateAxes::IntegrateAxes(SW_Builder& builderInput, double min, double max, double points)
{
    builder = builderInput;
    MinimumEnergy = min;
    MaximumEnergy = max;
    EnergyPoints = points;
}

int IntegrateAxes::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    TwoDimensionResolutionFunction tmp(builder, MinimumEnergy, MaximumEnergy, EnergyPoints);
    vector<double> val = tmp.getCut(x[0],ky,x[1]);

    for(int i=0;i!=EnergyPoints;i++)
     {
     retval[i] = val[i];
     //double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
     //cout << energy << " " << val[i] << " ";
     }
     //cout << endl;
    return 0;
}

int IntegrateAxes::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
{
    // Call non-static member function.
    return static_cast<IntegrateAxes*>(data)->calculateIntegrand(dim,x,fdim,retval);
}

std::vector<double> IntegrateAxes::getCut(double kxIn, double kyIn, double kzIn)
{
    kx = kxIn;
    ky = kyIn;
    kz = kzIn;
    
    vector<double> xmin(2);
    vector<double> xmax(2);
    
    xmin[0] = kx-0.2;
    xmax[0] = kx+0.2;
    xmin[1] = kz-0.2;
    xmax[1] = kz+0.2;
    
    double tol = 1.0e-4;
    
    vector<double> fval(EnergyPoints);
    vector<double> err(EnergyPoints);
    
    hcubature(EnergyPoints,IntegrateAxes::calc, this, 2, &xmin[0], &xmax[0], 0, tol, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    /*for(int i=0;i!=EnergyPoints;i++)
    {
        double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
        cout << energy << " " << fval[i] << " ";
    }
    cout << endl;
    */
    return fval;
    
}

/*NoResolutionFunction::NoResolutionFunction()
 {
 
 }
 
 std::vector<double> NoResolutionFunction::getCut(double kx, double ky, double kz)
 {
 }
 
 EnergyResolutionFunction::EnergyResolutionFunction()
 {
 }
 
 std::vector<double> EnergyResolutionFunction::getCut(double kx, double ky, double kz)
 {
 
 }
 
 TwoDimensionResolutionFunction::TwoDimensionResolutionFunction()
 {
 
 }*/

/*FourDimensionResolutionFunction::FourDimensionResolutionFunction()
 {
 
 }
 
 std::vector<double> FourDimensionResolutionFunction::getCut(double kx, double ky, double kz)
 {
 }
 
 IntegrateAxes::IntegrateAxes()
 {
 
 }*/