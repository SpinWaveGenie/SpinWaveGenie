//
//  TwoDimensionalGaussian.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/21/14.
//
//


#include <cmath>
#include <algorithm>
#include <functional>
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "External/cubature.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Plot/TwoDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/TwoDGaussian.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"

using namespace std;

namespace SpinWaveGenie
{
    TwoDimensionResolutionFunction::TwoDimensionResolutionFunction(TwoDimGaussian& info, SpinWave SW, Energies energiesIn)
    {
        a = info.a;
        b = info.b;
        c = info.c;
        info.direction.normalize();
        direction = info.direction;
        tolerance = info.tol;
        maximumEvaluations = 100000;
        energies = energiesIn;
        res.setSpinWave(SW);
        res.setEnergies(energies);
    }
    
    int TwoDimensionResolutionFunction::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *fval)
    {
        //cout << x[0] << " " << ky << " " << kz << " " << u << endl;
        unique_ptr<TwoDGaussian> resinfo(new TwoDGaussian());
        resinfo->setTolerance(0.01*tolerance);
        resinfo->setResolution(a,b,c);
        resinfo->setU(x[0]);
    
        res.setResolutionFunction(move(resinfo));
        vector<double> result = res.getCut(kx+x[0]*direction[0],ky+x[0]*direction[1],kz+x[0]*direction[2]);
        std::copy(result.begin(),result.end(),fval);
        return 0;
    }
    
    int TwoDimensionResolutionFunction::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
    {
    // Call non-static member function.
        return static_cast<TwoDimensionResolutionFunction*>(data)->calculateIntegrand(dim,x,fdim,retval);
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
    
    std::vector<double> TwoDimensionResolutionFunction::getCut(double kxIn, double kyIn, double kzIn)
    {
        size_t EnergyPoints = energies.size();
        kx = kxIn;
        ky = kyIn;
        kz = kzIn;
    
        double xmin, xmax;
        double d = -log(tolerance);
        xmax = sqrt(c*d/(a*c-b*b));
        xmin = -xmax;
    
        vector<double> fval(EnergyPoints);
        vector<double> err(EnergyPoints);
    
        double tmp = sqrt((M_PI*M_PI)/(a*c-b*b));
        hcubature(EnergyPoints,TwoDimensionResolutionFunction::calc, this, 1, &xmin, &xmax, maximumEvaluations, tolerance/tmp, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
        std::for_each(fval.begin(), fval.end(), DivideValue(tmp));
        return fval;
    }
    
    const Cell& TwoDimensionResolutionFunction::getCell() const
    {
        return res.getCell();
    }
    
    const Energies& TwoDimensionResolutionFunction::getEnergies()
    {
        return energies;
    }
    
    void TwoDimensionResolutionFunction::setTolerance(double toleranceIn, int maxEvals)
    {
        tolerance = toleranceIn;
        maximumEvaluations = maxEvals;
    }
    
    void TwoDimensionResolutionFunction::setEnergies(Energies energiesIn)
    {
        energies = energiesIn;
    }
    
    std::unique_ptr<SpinWavePlot> TwoDimensionResolutionFunction::clone()
    {
        return unique_ptr<SpinWavePlot>(new TwoDimensionResolutionFunction(*this));
    }
    
}
