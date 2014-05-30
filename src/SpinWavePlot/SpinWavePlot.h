#ifndef __SpinWavePlot__
#define __SpinWavePlot__


#include <iostream>
#include <string>
#include <vector>
#include "Cell/Cell.h"
#include "Genie/SpinWave.h"
#include "OneDimensionalGaussian.h"
#include "Containers/Energies.h"

/* Abstract base class */
class SpinWavePlot
{
public:
    virtual std::unique_ptr<SpinWavePlot> clone() = 0;
    virtual const Cell& getCell() const = 0;
    virtual const Energies& getEnergies() = 0;
    virtual void setEnergies(Energies energies) = 0;
    virtual std::vector<double> getCut(double kx, double ky, double kz) = 0; // returns constant-Q cut
    virtual ~SpinWavePlot(){};
};

#endif /* defined(__SpinWavePlot__) */
