#ifndef __SpinWavePlot__
#define __SpinWavePlot__


#include <iostream>
#include <string>
#include <vector>
#include "Cell/Cell.h"
#include "Genie/SpinWave.h"
#include "OneDimensionalGaussian.h"

/* Abstract base class */
class SpinWavePlot
{
public:
    virtual std::unique_ptr<SpinWavePlot> clone() = 0;
    virtual const Cell& getCell() const = 0;
    virtual double getMinimumEnergy() const = 0;
    virtual void setMinimumEnergy(double energy) = 0;
    virtual double getMaximumEnergy() const = 0;
    virtual void setMaximumEnergy(double energy) = 0;
    virtual std::size_t getNumberPoints() const = 0;
    virtual void setNumberPoints(std::size_t points) = 0;
    virtual std::vector<double> getCut(double kx, double ky, double kz) = 0; // returns constant-Q cut
    virtual ~SpinWavePlot(){};
};

struct TwoDimGaussian
{
    double a,b,c;
    unsigned direction;
    double tol;
    SpinWave SW;
};

class TwoDimensionResolutionFunction : public SpinWavePlot{
public:
    TwoDimensionResolutionFunction(){};
    TwoDimensionResolutionFunction(TwoDimGaussian& info, double min, double max, double points);
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
    ~TwoDimensionResolutionFunction(){};
private:
    double MinimumEnergy,MaximumEnergy,tol;
    double a,b,c;
    double kx,ky,kz;
    unsigned EnergyPoints,direction;
    SpinWave SW;
};




/*
 class FourDimensionResolutionFunction : SpinWavePlot {
 FourDimensionResolutionFunction();
 std::vector<double> getCut(double kx, double ky, double kz);
 ~FourDimensionResolutionFunction() = 0;
 };
 */

#endif /* defined(__SpinWavePlot__) */
