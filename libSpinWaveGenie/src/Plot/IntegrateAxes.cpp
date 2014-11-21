//
//  IntegrateAxes.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/7/14.
//
//
#include "SpinWaveGenie/Plot/IntegrateAxes.h"
#include "cuba.h"

using namespace std;

namespace SpinWaveGenie
{

IntegrateAxes::IntegrateAxes(const IntegrateAxes& other)
{
    resolutionFunction = move(other.resolutionFunction->clone());
    this->maximumEvaluations = other.maximumEvaluations;
    this->tolerance = other.tolerance;
    this->integrationDirections = other.integrationDirections;
}

IntegrateAxes::IntegrateAxes(unique_ptr<SpinWavePlot> resFunction, HKLDirections directions, double tol, int maxEvals)
{
    setenv("CUBACORES","0",1);
    this->tolerance = tol;
    this->maximumEvaluations = maxEvals;
    this->integrationDirections = directions;
    this->resolutionFunction = move(resFunction);
}

int IntegrateAxes::calculateIntegrand(const int* ndim, const double *x, const int* ncomp, double *retval)
{
    int dim = *ndim;
    double tmpx=kx,tmpy=ky,tmpz=kz;
    //cout << dim << " dimensions" << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    
    for(unsigned i = 0; i!=dim; i++)
    {
        double xx = xmin[i]+x[i]*(xmax[i]-xmin[i]);
        //cout << i << "\t" << xx << endl;
        tmpx += xx*integrationDirections[i].v0;
        tmpy += xx*integrationDirections[i].v1;
        tmpz += xx*integrationDirections[i].v2;
    }
    //cout << endl;
    
    //cout << "** " << x[0] << endl;//<< " " << x[1] << " " << x[2] << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    vector<double> val = resolutionFunction->getCut(tmpx,tmpy,tmpz);
    
    std::copy(val.begin(),val.end(),retval);
    
    return 0;
}

int IntegrateAxes::calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
{

    return static_cast<IntegrateAxes*>(userdata)->calculateIntegrand(ndim,xx,ncomp,ff);
}

std::vector<double> IntegrateAxes::getCut(double kxIn, double kyIn, double kzIn)
{
    kx = kxIn;
    ky = kyIn;
    kz = kzIn;
    
    int dim = integrationDirections.size();
    
    xmin.clear();
    xmin.reserve(dim);
    xmax.clear();
    xmax.reserve(dim);
    
    for (auto it = integrationDirections.begin(); it != integrationDirections.end(); it++)
    {
        //cout << -1.0*it->delta << " " << it->delta << endl;
        xmin.push_back(-1.0*it->delta);
        xmax.push_back(it->delta);
    }
    
    //cout << "** " << kx << " " << ky << " " << kz << endl;
    //cout << xmin[0] << " " << endl; //<< xmin[1] << " " << xmin[2] << endl;
    //cout << xmax[0] << " " << endl;//<< xmax[1] << " " << xmax[2] << endl;
    size_t energyPoints = resolutionFunction->getEnergies().size();    
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    vector<double> prob(energyPoints);
    
    int nregions, // the actual number of subregions needed
    neval, // the actual number of integrand evaluations needed.
    fail; // error flag: 0, the desired accuracy was reached. −1, dimension out of range. > 0, the accuracy goal was not met within the allowed maximum number of integrand evaluations.
    Cuhre(dim,energyPoints, IntegrateAxes::calc, this, 1, 1000.0, tolerance, 0, 0, 50000,0,NULL, NULL,&nregions, &neval, &fail, &fval[0],&err[0], &prob[0]);
    return fval;
}

const Cell& IntegrateAxes::getCell() const
{
    return resolutionFunction->getCell();
}

const Energies& IntegrateAxes::getEnergies()
{
    return resolutionFunction->getEnergies();
}

void IntegrateAxes::setEnergies(Energies energiesIn)
{
    resolutionFunction->setEnergies(energiesIn);
}

std::unique_ptr<SpinWavePlot> IntegrateAxes::clone()
{
    return unique_ptr<SpinWavePlot>(new IntegrateAxes(*this));
}
    
}
