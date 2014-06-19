#include <cmath>
#include "OneDimensionalGaussian.h"

using namespace std;

namespace SpinWaveGenie
{

void OneDimensionalGaussian::setFWHM(double InFWHM)
{
    FWHM = InFWHM;
}

void OneDimensionalGaussian::setTolerance(double InTolerance)
{
    Tolerance = InTolerance;
}

double OneDimensionalGaussian::getMinimumEnergy()
{
    return -1.0*getMaximumEnergy();
}

double OneDimensionalGaussian::getMaximumEnergy()
{
    double factor = 2.0*sqrt(log(2.0))/(FWHM*sqrt(M_PI));
    double constant = log(Tolerance/factor);
    double ma = getExponentialFactor();
    return sqrt(constant/ma);
}

double OneDimensionalGaussian::getFunction(double frequency, double energy)
{
    //cout << "frequency: " << frequency << endl;
    //cout << "energy: " << energy << endl;
    double ma = getExponentialFactor();
    double factor = 2.0*sqrt(log(2.0))/(FWHM*sqrt(M_PI));
    double result = factor*exp(ma*pow(frequency-energy,2));
    //cout << "result:" << ma << " " << factor << " " << result << endl;
    return result;
}

double OneDimensionalGaussian::getExponentialFactor()
{
    return -4.0*log(2.0)/pow(FWHM,2);
}

unique_ptr<OneDimensionalShapes> OneDimensionalGaussian::clone()
{
    return unique_ptr<OneDimensionalShapes>(new OneDimensionalGaussian(*this));
}
}