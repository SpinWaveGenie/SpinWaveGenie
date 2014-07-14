//
//  OneDimensionalFactory.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/10/14.
//
//

#ifndef __OneDimensionalFactory__
#define __OneDimensionalFactory__

#include <iostream>
#include <memory>
#include "SpinWaveGenie/SpinWavePlot/OneDimensionalShapes.h"

namespace SpinWaveGenie
{

class OneDimensionalFactory
{
public:
    std::unique_ptr<OneDimensionalShapes> getGaussian(double fwhm, double tol);
    std::unique_ptr<OneDimensionalShapes> getLorentzian(double fwhm, double tol);
    std::unique_ptr<OneDimensionalShapes> getPseudoVoigt(double eta, double fwhm, double tol);
private:
};
}
#endif /* defined(__OneDimensionalFactory__) */
