//
//  OneDimensionalPseudoVoigt.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/6/14.
//
//

#ifndef __OneDimensionalPseudoVoigt__
#define __OneDimensionalPseudoVoigt__

#include <iostream>
#include <memory>
#include "OneDimensionalShapes.h"
#include "OneDimensionalLorentzian.h"
#include "OneDimensionalGaussian.h"

class OneDimensionalPseudoVoigt: public OneDimensionalShapes
{
public:
    OneDimensionalPseudoVoigt();
    OneDimensionalPseudoVoigt(const OneDimensionalPseudoVoigt& other) : Lorentzian( other.Lorentzian->clone() ) , Gaussian( other.Gaussian->clone() ) {};
    void setEta(double InEta);
    void setFWHM(double InFWHM);
    void setTolerance(double InTolerance);
    double getMinimumEnergy();
    double getMaximumEnergy();
    double getFunction(double frequency, double energy);
    std::unique_ptr<OneDimensionalShapes> clone();
    ~OneDimensionalPseudoVoigt(){};
private:
    double eta;
    std::unique_ptr<OneDimensionalShapes> Lorentzian,Gaussian;
};

#endif /* defined(__OneDimensionalPseudoVoigt__) */
