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

IntegrateAxes::IntegrateAxes(unique_ptr<SpinWavePlot> resFunction, double tolerance)
{
    this->tolerance = tolerance;
    this->resolutionFunction = move(resFunction);
    this->minimumEnergy = resolutionFunction->getMinimumEnergy();
    this->maximumEnergy = resolutionFunction->getMaximumEnergy();
    this->energyPoints = resolutionFunction->getNumberPoints();
}

void IntegrateAxes::addDirection(double v0, double v1, double v2, double delta)
{
    Axis direction;
    direction.v0 = v0;
    direction.v1 = v1;
    direction.v2 = v2;
    direction.delta = delta;
    
    integrationDirections.push_back(direction);
}

void IntegrateAxes::addDirection(int direction, double delta)
{
    double v0 = 0.0, v1=0.0, v2=0.0;
    switch(direction)
    {
        case 0:
            v0 = 1.0;
            break;
        case 1:
            v1 = 1.0;
            break;
        case 2:
            v2 = 1.0;
            break;
    }
    
    this->addDirection(v0, v1, v2, delta);
}

int IntegrateAxes::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
{
    double tmpx=kx,tmpy=ky,tmpz=kz;
    //cout << dim << " dimensions" << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    switch(dim)
    {
        case 1:
        {
            tmpx += x[0]*integrationDirections[0].v0;
            tmpy += x[0]*integrationDirections[0].v1;
            tmpz += x[0]*integrationDirections[0].v2;
            break;
        }
        case 2:
        {
            tmpx += x[0]*integrationDirections[0].v0+x[1]*integrationDirections[1].v0;
            tmpy += x[0]*integrationDirections[0].v1+x[1]*integrationDirections[1].v1;
            tmpz += x[0]*integrationDirections[0].v2+x[1]*integrationDirections[1].v2;
            break;
        }
        case 3:
        {
            tmpx += x[0]*integrationDirections[0].v0+x[1]*integrationDirections[1].v0+x[2]*integrationDirections[2].v0;
            tmpy += x[0]*integrationDirections[0].v1+x[1]*integrationDirections[1].v1+x[2]*integrationDirections[2].v1;
            tmpz += x[0]*integrationDirections[0].v2+x[1]*integrationDirections[1].v2+x[2]*integrationDirections[2].v2;
            break;
        }
    }
    //cout << "** " << x[0] << endl;//<< " " << x[1] << " " << x[2] << endl;
    //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
    
    vector<double> val = resolutionFunction->getCut(tmpx,tmpy,tmpz);
    
    std::copy(val.begin(),val.end(),retval);

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
    
    std::vector<double> xmin,xmax;
    
    
    for (auto it = integrationDirections.begin(); it != integrationDirections.end(); it++)
    {
        xmin.push_back(-1.0*it->delta);
        xmax.push_back(it->delta);
    }
    
    Matrix3 lattice = resolutionFunction->getCell().getReciprocalVectors();
    int dim = integrationDirections.size();

    switch(dim)
    {
        case 1:
        {
            Vector3 direction = integrationDirections[0].v0*lattice.row(0) + integrationDirections[0].v1*lattice.row(1) + integrationDirections[0].v2*lattice.row(2);
            direction *= 2.0*integrationDirections[0].delta;
            volume = direction.norm();
            break;
        }
        case 2:
        {
            Vector3 firstDirection = integrationDirections[0].v0*lattice.row(0) + integrationDirections[0].v1*lattice.row(1) + integrationDirections[0].v2*lattice.row(2);
            firstDirection *= 2.0*integrationDirections[0].delta;
            Vector3 secondDirection = integrationDirections[1].v0*lattice.row(0) + integrationDirections[1].v1*lattice.row(1) + integrationDirections[1].v2*lattice.row(2);
            secondDirection *= 2.0*integrationDirections[1].delta;
            volume = firstDirection.cross(secondDirection).norm();
            break;
        }
        case 3:
        {
            Vector3 firstDirection = integrationDirections[0].v0*lattice.row(0) + integrationDirections[0].v1*lattice.row(1) + integrationDirections[0].v2*lattice.row(2);
            firstDirection *= 2.0*integrationDirections[0].delta;
            Vector3 secondDirection = integrationDirections[1].v0*lattice.row(0) + integrationDirections[1].v1*lattice.row(1) + integrationDirections[1].v2*lattice.row(2);
            secondDirection *= 2.0*integrationDirections[1].delta;
            Vector3 thirdDirection = integrationDirections[2].v0*lattice.row(0) + integrationDirections[2].v1*lattice.row(1) + integrationDirections[2].v2*lattice.row(2);
            thirdDirection *= 2.0*integrationDirections[2].delta;

            volume = std::abs(firstDirection.dot(secondDirection.cross(thirdDirection)));
            break;
        }
    }
    
    //cout << "** " << kx << " " << ky << " " << kz << endl;
    //cout << xmin[0] << " " << endl; //<< xmin[1] << " " << xmin[2] << endl;
    //cout << xmax[0] << " " << endl;//<< xmax[1] << " " << xmax[2] << endl;
    
    vector<double> fval(energyPoints);
    vector<double> err(energyPoints);
    
    hcubature(energyPoints,IntegrateAxes::calc, this, dim, &xmin[0], &xmax[0], 0, tolerance, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
    
    //cout << "volume = " << volume << endl;

    std::for_each(fval.begin(), fval.end(), DivideValue(volume));

    return fval;
}