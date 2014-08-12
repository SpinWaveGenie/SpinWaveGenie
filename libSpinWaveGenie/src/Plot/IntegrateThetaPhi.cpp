//
//  IntegrateThetaPhi.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#include "SpinWaveGenie/Plot/IntegrateThetaPhi.h"
#include <vector>
#include "External/cubature.h"


using namespace std;

namespace SpinWaveGenie
{

IntegrateThetaPhi::IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tolerance)
{
    tol = tolerance;
    resolutionFunction = move(object);
}

IntegrateThetaPhi::IntegrateThetaPhi(const IntegrateThetaPhi& other)
{
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

const Energies& IntegrateThetaPhi::getEnergies()
{
    return resolutionFunction->getEnergies();
}

void IntegrateThetaPhi::setEnergies(Energies energiesIn)
{
    resolutionFunction->setEnergies(energiesIn);
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
    //cout << tmp.norm() << endl;
    //cout << k.transpose() << endl;
    
    vector<double> val = resolutionFunction->getCut(k[0],k[1],k[2]);
    
    //for (auto it = values.begin(); it!=values.end(); it++)
    //{
    //    cout << (*it) << endl;
    //}
    
    //cout << MinimumEnergy << " " << MaximumEnergy << " " << EnergyPoints << endl;
    double factor = sin(theta)/(4.0*M_PI);
    size_t energyPoints = resolutionFunction->getEnergies().size();
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
    size_t energyPoints = resolutionFunction->getEnergies().size();
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    
    hcubature(energyPoints,IntegrateThetaPhi::calc, this, dim, xmin.data(), xmax.data(), 100000, tol, 0, ERROR_INDIVIDUAL, fval.data(), err.data());
    
    /*for(int i=0;i!=EnergyPoints;i++)
     {
     double energy = MinimumEnergy + (MaximumEnergy-MinimumEnergy)*(double)i/(double)(EnergyPoints-1);
     cout << energy << " " << fval[i] << " ";
     }
     cout << endl;
     */
    return fval;
}
}
