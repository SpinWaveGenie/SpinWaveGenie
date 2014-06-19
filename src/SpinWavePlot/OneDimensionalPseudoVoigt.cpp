//
//  OneDimensionalPseudoVoigt.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/6/14.
//
//

#include <cmath>
#include "OneDimensionalPseudoVoigt.h"

using namespace std;

namespace SpinWaveGenie
{

OneDimensionalPseudoVoigt::OneDimensionalPseudoVoigt()
{
    Gaussian.reset(new OneDimensionalGaussian);
    Lorentzian.reset(new OneDimensionalLorentzian);
}

void OneDimensionalPseudoVoigt::setEta(double InEta)
{
    eta = InEta;
}

void OneDimensionalPseudoVoigt::setFWHM(double InFWHM)
{
    Gaussian->setFWHM(InFWHM);
    Lorentzian->setFWHM(InFWHM);
}

void OneDimensionalPseudoVoigt::setTolerance(double InTolerance)
{
    Gaussian->setTolerance(InTolerance);
    Lorentzian->setTolerance(InTolerance);
}

double OneDimensionalPseudoVoigt::getMinimumEnergy()
{
    return -1.0*getMaximumEnergy();
}

double OneDimensionalPseudoVoigt::getMaximumEnergy()
{
    return Lorentzian->getMaximumEnergy();
}

double OneDimensionalPseudoVoigt::getFunction(double frequency, double energy)
{
    double result = eta*Lorentzian->getFunction(frequency,energy) + (1.0-eta)*Gaussian->getFunction(frequency,energy);
    return result;
}

unique_ptr<OneDimensionalShapes> OneDimensionalPseudoVoigt::clone()
{
    return unique_ptr<OneDimensionalShapes>(new OneDimensionalPseudoVoigt(*this));
}

}