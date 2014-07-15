//
//  TwoDGaussian.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/23/14.
//
//

#ifndef __spin_wave_genie__TwoDGaussian__
#define __spin_wave_genie__TwoDGaussian__

#include <iostream>
#include <memory>
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"

namespace SpinWaveGenie
{
    
    class TwoDGaussian: public OneDimensionalShapes
    {
    public:
        //void setFWHM(double InFWHM){};
        void setResolution( double aIn, double bIn, double cIn);
        void setU( double uIn);
        void setTolerance(double InTolerance);
        double getMinimumEnergy();
        double getMaximumEnergy();
        double getFunction(double frequency, double energy);
        std::unique_ptr<OneDimensionalShapes> clone();
        ~TwoDGaussian(){};
    private:
        double getExponentialFactor();
        double FWHM,Tolerance;
        double a,b,c;
        double u;
    };
}

#endif /* defined(__spin_wave_genie__TwoDGaussian__) */
