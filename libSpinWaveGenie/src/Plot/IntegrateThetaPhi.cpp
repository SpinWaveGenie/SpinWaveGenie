//
//  IntegrateThetaPhi.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#include "SpinWaveGenie/Plot/IntegrateThetaPhi.h"
#include <vector>
#include "cuba.h"


using namespace std;

namespace SpinWaveGenie
{

IntegrateThetaPhi::IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tol, int maxEvals)
{
    setenv("CUBACORES","0",1);
    tolerance = tol;
    maximumEvaluations = maxEvals;
    resolutionFunction = move(object);
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


int IntegrateThetaPhi::calculateIntegrand(const int* dim, const double *x,const int* fdim, double *retval)
{
    Vector3 tmp,k;
    double theta = x[0]*M_PI;
    double phi = x[1]*2.0*M_PI;
    
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
    double factor = sin(theta)*M_PI_2;
    size_t energyPoints = resolutionFunction->getEnergies().size();
    for(size_t i=0;i!=energyPoints;i++)
    {
        retval[i] = factor*val[i];
    }
    //cout << endl;
    return 0;
}

int IntegrateThetaPhi::calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{
    // Call non-static member function.
    return static_cast<IntegrateThetaPhi*>(userdata)->calculateIntegrand(ndim,xx,ncomp,ff);
}

std::vector<double> IntegrateThetaPhi::getCut(double kx,double ky, double kz)
{
    //std::vector<double> xmin = {0.0,0.0};
    //std::vector<double> xmax = {M_PI,2.0*M_PI};
    int dim = 2;
    r = std::abs(kz);
    
    //cout << "dispAng = " << r << endl;
    size_t energyPoints = resolutionFunction->getEnergies().size();
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    vector<double> prob(energyPoints);
    
    int nregions, // the actual number of subregions needed
    neval, // the actual number of integrand evaluations needed.
    fail; // error flag: 0, the desired accuracy was reached. −1, dimension out of range. > 0, the accuracy goal was not met within the allowed maximum number of integrand evaluations.
    Cuhre(dim,energyPoints, IntegrateThetaPhi::calc, this, 1, 1000.0, tolerance, 0, 0, 50000,0,NULL, NULL,&nregions, &neval, &fail, &fval[0],&err[0], &prob[0]);

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
