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
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Plot/TwoDimensionalGaussian.h"
#include "SpinWaveGenie/Plot/TwoDGaussian.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"
#include "AdaptiveSimpson.h"

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
        maximumEvaluations = 100;
        energies = energiesIn;
        res.setSpinWave(SW);
        res.setEnergies(energies);
    }
    
    std::vector<double> TwoDimensionResolutionFunction::calculateIntegrand(std::deque<double>& x)
    {
        unique_ptr<TwoDGaussian> resinfo(new TwoDGaussian());
        resinfo->setTolerance(0.1*tolerance);
        resinfo->setResolution(a,b,c);
        resinfo->setU(x[0]);
    
        res.setResolutionFunction(move(resinfo));
        return res.getCut(kx+x[0]*direction[0],ky+x[0]*direction[1],kz+x[0]*direction[2]);
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
        kx = kxIn;
        ky = kyIn;
        kz = kzIn;
    
        std::vector<double> xmin(1), xmax(1);
        double d = -log(tolerance);
        xmax[0] = sqrt(c*d/(a*c-b*b));
        xmin[0] = -1.0*xmax[0];

        double tmp = sqrt((M_PI*M_PI)/(a*c-b*b));

        std::function< std::vector<double>(std::deque<double>& x)> funct = std::bind<std::vector<double> >(&TwoDimensionResolutionFunction::calculateIntegrand,this,std::placeholders::_1);
        AdaptiveSimpson test;
        test.setFunction(funct);
        test.setInterval(xmin,xmax);
        test.setPrecision(tolerance);
        test.setMaximumRecursionDepth(maximumEvaluations);
        std::vector<double> result = test.integrate();
        std::for_each(result.begin(), result.end(), DivideValue(tmp));
        return result;
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
