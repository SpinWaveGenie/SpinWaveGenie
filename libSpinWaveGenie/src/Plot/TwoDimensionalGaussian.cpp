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
#include "cuba.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Plot/TwoDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/TwoDGaussian.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"

using namespace std;

namespace SpinWaveGenie
{
    TwoDimensionResolutionFunction::TwoDimensionResolutionFunction(TwoDimGaussian& info, SpinWave SW, Energies energiesIn)
    {
        setenv("CUBACORES","0",1);
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
    
    int TwoDimensionResolutionFunction::calculateIntegrand(const int* ndim, const double *x, const int* ncomp, double *retval)
    {
        //int dim = *ndim;
        //int fdim = *ncomp;
        //cout << x[0] << " " << ky << " " << kz << " " << u << endl;
        double scaled_u = x[0]*2.0*xmax;
        
        unique_ptr<TwoDGaussian> resinfo(new TwoDGaussian());
        resinfo->setTolerance(0.01*tolerance);
        resinfo->setResolution(a,b,c);
        resinfo->setU(scaled_u);
    
        res.setResolutionFunction(move(resinfo));
        vector<double> result = res.getCut(kx+scaled_u*direction[0],ky+scaled_u*direction[1],kz+scaled_u*direction[2]);
        std::copy(result.begin(),result.end(),retval);
        return 0;
    }
    
    int TwoDimensionResolutionFunction::calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
    {
        // Call non-static member function.
        return static_cast<TwoDimensionResolutionFunction*>(userdata)->calculateIntegrand(ndim,xx,ncomp,ff);

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
    
        double d = -log(tolerance);
        xmax = sqrt(c*d/(a*c-b*b));
    
        vector<double> fval(EnergyPoints);
        vector<double> err(EnergyPoints);
        vector<double> prob(EnergyPoints);
        
        int nregions, neval, fail;
        Cuhre(1,EnergyPoints, TwoDimensionResolutionFunction::calc, this, 1, 1000.0, tolerance*2.0*xmax, 0, 0, 50000,0,NULL, NULL,&nregions, &neval, &fail, &fval[0],&err[0], &prob[0]);
        double tmp = sqrt((M_PI*M_PI)/(a*c-b*b))/(2.0*xmax);
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
