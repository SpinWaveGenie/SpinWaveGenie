//
//  IntegrateAxes.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/7/14.
//
//
#include "IntegrateAxes.h"
#include "External/cubature.h"

using namespace std;

IntegrateAxes::IntegrateAxes(axes_info info, unique_ptr<SpinWavePlot> resFunction)
{
    this->info = info;
    
    ResolutionFunction = move(resFunction);
    MinimumEnergy = ResolutionFunction->getMinimumEnergy();
    MaximumEnergy = ResolutionFunction->getMaximumEnergy();
    EnergyPoints = ResolutionFunction->getNumberPoints();
}

int IntegrateAxes::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    double tmpx=kx,tmpy=ky,tmpz=kz;
    
    size_t placeholder = 0;
    if(info.x)
    {
        tmpx += x[placeholder];
        placeholder++;
    }
    
    if(info.y)
    {
        tmpx += 0.5*x[placeholder];
        tmpy -= x[placeholder];
        placeholder++;
    }

    if(info.z)
    {
        tmpz += x[placeholder];
        placeholder++;
    }
    
    //cout << "** " << x[0] << " " << x[1] << " " << x[2] << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    vector<double> val = ResolutionFunction->getCut(tmpx,tmpy,tmpz);
    
    for(int i=0;i!=EnergyPoints;i++)
    {
        retval[i] = val[i]/volume;
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
    
    std::vector<double> xmin,xmax;
    
    Matrix3 lattice = ResolutionFunction->getCell().getReciprocalVectors();
    int dim = 0;
    volume = 1.0;
    if (info.x)
    {
        dim++;
        xmin.push_back(-1.0*info.dx);
        xmax.push_back(info.dx);
        double tmp = lattice.row(0).norm();
        volume *= 2.0*tmp*info.dx;
    }
    if (info.y)
    {
        dim++;
        xmin.push_back(-1.0*info.dy);
        xmax.push_back(info.dy);
        double tmp1 = lattice.row(0).norm();
        double tmp2 = lattice.row(1).norm();
        volume *= sqrt(tmp1*tmp1+4.0*tmp2*tmp2)*info.dy;
    
    }
    if (info.z)
    {
        dim++;
        xmin.push_back(-1.0*info.dz);
        xmax.push_back(info.dz);
        double tmp = lattice.row(2).norm();
        volume *= 2.0*tmp*info.dz;
    }
    
    //cout << "** " << kx << " " << ky << " " << kz << endl;
    //cout << xmin[0] << " " << xmin[1] << " " << xmin[2] << endl;
    //cout << xmax[0] << " " << xmax[1] << " " << xmax[2] << endl;
    
    vector<double> fval(EnergyPoints);
    vector<double> err(EnergyPoints);
    
    hcubature(EnergyPoints,IntegrateAxes::calc, this, dim, &xmin[0], &xmax[0], 0, info.tol, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    /*for(int i=0;i!=EnergyPoints;i++)
     {
     double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
     cout << energy << " " << fval[i] << " ";
     }
     cout << endl;
     */
    return fval;
}