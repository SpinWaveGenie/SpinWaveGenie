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


IntegrateThetaPhi::IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tolerance)
{
    minimumEnergy = object->getMinimumEnergy();
    maximumEnergy = object->getMaximumEnergy();
    energyPoints = object->getNumberPoints();
    tol = tolerance;
    resolutionFunction = move(object);
}

IntegrateThetaPhi::IntegrateThetaPhi(const IntegrateThetaPhi& other)
{
    minimumEnergy = other.minimumEnergy;
    maximumEnergy = other.maximumEnergy;
    energyPoints = other.energyPoints;
    tol = other.tol;
    resolutionFunction = move(other.resolutionFunction->clone());
}

std::unique_ptr<SpinWavePlot> IntegrateThetaPhi::clone()
{
    return unique_ptr<SpinWavePlot>(new IntegrateThetaPhi(*this));
}

const Cell& IntegrateThetaPhi::getCell() const
{
    return resolutionFunction->getCell();
}

double IntegrateThetaPhi::getMinimumEnergy() const
{
    return minimumEnergy;
}

void IntegrateThetaPhi::setMinimumEnergy(double energy)
{
    resolutionFunction->setMinimumEnergy(energy);
    this->minimumEnergy = energy;
}

double IntegrateThetaPhi::getMaximumEnergy() const
{
    return maximumEnergy;
}

void IntegrateThetaPhi::setMaximumEnergy(double energy)
{
    resolutionFunction->setMaximumEnergy(energy);
    this->maximumEnergy = energy;
}

std::size_t IntegrateThetaPhi::getNumberPoints() const
{
    return energyPoints;
}

void IntegrateThetaPhi::setNumberPoints(std::size_t points)
{
    resolutionFunction->setNumberPoints(points);
    this->energyPoints = points;
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
    
    Matrix3 basisVectors = resolutionFunction->getCell().getBasisVectors();
    
    k = tmp.transpose()*basisVectors/(2.0*M_PI);
    //cout << tmp.transose() << endl;
    //cout << k.transpose() << endl;
    
    vector<double> val = resolutionFunction->getCut(k[0],k[1],k[2]);
    
    //for (auto it = values.begin(); it!=values.end(); it++)
    //{
    //    cout << (*it) << endl;
    //}
    
    //cout << MinimumEnergy << " " << MaximumEnergy << " " << EnergyPoints << endl;
    double factor = sin(theta)/(4.0*M_PI);
    for(int i=0;i!=energyPoints;i++)
    {
        retval[i] = factor*val[i];
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
    r = std::abs(kz);
    
    //cout << "dispAng = " << r << endl;
    
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    
    hcubature(energyPoints,IntegrateThetaPhi::calc, this, dim, &xmin[0], &xmax[0], 0, tol, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    /*for(int i=0;i!=EnergyPoints;i++)
     {
     double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
     cout << energy << " " << fval[i] << " ";
     }
     cout << endl;
     */
    return fval;
}
