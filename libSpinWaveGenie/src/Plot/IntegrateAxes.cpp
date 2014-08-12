//
//  IntegrateAxes.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/7/14.
//
//
#include "SpinWaveGenie/Plot/IntegrateAxes.h"
#include "External/cubature.h"

using namespace std;

namespace SpinWaveGenie
{

IntegrateAxes::IntegrateAxes(const IntegrateAxes& other)
{
    resolutionFunction = move(other.resolutionFunction->clone());
    this->tolerance = other.tolerance;
    this->integrationDirections = other.integrationDirections;
}

IntegrateAxes::IntegrateAxes(unique_ptr<SpinWavePlot> resFunction, HKLDirections directions, double tolerance)
{
    this->tolerance = tolerance;
    this->integrationDirections = directions;
    this->resolutionFunction = move(resFunction);
}

int IntegrateAxes::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    double tmpx=kx,tmpy=ky,tmpz=kz;
    //cout << dim << " dimensions" << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    for(unsigned i = 0; i!=dim; i++)
    {
        tmpx += x[i]*integrationDirections[i].v0;
        tmpy += x[i]*integrationDirections[i].v1;
        tmpz += x[i]*integrationDirections[i].v2;
    }
    
    
    //cout << "** " << x[0] << endl;//<< " " << x[1] << " " << x[2] << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    vector<double> val = resolutionFunction->getCut(tmpx,tmpy,tmpz);
    
    std::copy(val.begin(),val.end(),retval);

    for (int i=0;i<val.size();i++)
    {
        if (retval[i] < 0.0)
            cout << "0.0 > fval[i] = " << retval[i] << endl;
    }
    
    return 0;
}

int IntegrateAxes::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
{
    // Call non-static member function.
    return static_cast<IntegrateAxes*>(data)->calculateIntegrand(dim,x,fdim,retval);
}


struct DivideValue
{
    double value;
    DivideValue(double v)
    {
        value = 1.0/v;
    }
    void operator()(double &elem) const
    {
        elem *= value;
    }
};


std::vector<double> IntegrateAxes::getCut(double kxIn, double kyIn, double kzIn)
{
    kx = kxIn;
    ky = kyIn;
    kz = kzIn;
    
    int dim = integrationDirections.size();
    
    std::vector<double> xmin,xmax;
    xmin.reserve(dim);
    xmax.reserve(dim);
    
    for (auto it = integrationDirections.begin(); it != integrationDirections.end(); it++)
    {
        //cout << -1.0*it->delta << " " << it->delta << endl;
        xmin.push_back(-1.0*it->delta);
        xmax.push_back(it->delta);
    }
    
    volume = 1.0;
    for (auto it = integrationDirections.begin(); it!= integrationDirections.end(); ++it)
    {
        volume *= 2.0*it->delta;
    }
            
    //cout << volume << endl;
    
    //cout << "** " << kx << " " << ky << " " << kz << endl;
    //cout << xmin[0] << " " << endl; //<< xmin[1] << " " << xmin[2] << endl;
    //cout << xmax[0] << " " << endl;//<< xmax[1] << " " << xmax[2] << endl;
    size_t energyPoints = resolutionFunction->getEnergies().size();    
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    
    hcubature(energyPoints,IntegrateAxes::calc, this, dim, &xmin[0], &xmax[0], 100000, tolerance/volume, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    //cout << "volume = " << volume << endl;

    for (int i=0;i<fval.size();i++)
    {
        if (fval[i] < 0.0)
            cout << "before: 0.0 > fval[i] = " << fval[i] << endl;
    }
    
    std::for_each(fval.begin(), fval.end(), DivideValue(volume));

    /*for (int i=0;i<fval.size();i++)
    {
        if (fval[i] < 0.0)
            cout << "after: 0.0 > fval[i] = " << fval[i] << endl;
    }*/
    
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
