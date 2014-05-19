//
//  IntegrateThetaPhi.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 12/5/13.
//
//

#ifndef __spin_wave_genie__IntegrateThetaPhi__
#define __spin_wave_genie__IntegrateThetaPhi__

#include <iostream>
#include "SpinWavePlot.h"

class IntegrateThetaPhi : public SpinWavePlot {
public:
    IntegrateThetaPhi(std::unique_ptr<SpinWavePlot> object, double tolerance);
    IntegrateThetaPhi(const IntegrateThetaPhi& other);
    std::unique_ptr<SpinWavePlot> clone();
    const Cell& getCell() const;
    double getMinimumEnergy() const;
    void setMinimumEnergy(double energy);
    double getMaximumEnergy() const;
    void setMaximumEnergy(double energy);
    std::size_t getNumberPoints() const;
    void setNumberPoints(std::size_t points);
    std::vector<double> getCut(double kx,double ky, double kz);
    ~IntegrateThetaPhi(){};
private:
    int calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval);
    static int calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval);
    double minimumEnergy,maximumEnergy;
    double r;
    unsigned energyPoints;
    double tol;
    std::unique_ptr<SpinWavePlot> resolutionFunction;
};

#endif /* defined(__spin_wave_genie__IntegrateThetaPhi__) */
