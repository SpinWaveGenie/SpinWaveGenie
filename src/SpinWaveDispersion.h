#ifndef __SpinWaveDispersion__
#define __SpinWaveDispersion__

#include <iostream>
#include "Containers/ThreeVectors.h"
#include "Genie/SpinWave.h"

class SpinWaveDispersion
{
public:
    enum class Options {PrintPosition,PrintFrequency,PrintIntensity};
    void setOptions(Options PrintOptions, bool Value);
    void setFilename(std::string name);
    void setPoints(ThreeVectors<double> pos);
    void setGenie(SpinWave SW);
    void save();
private:
    SpinWave Genie;
    std::string Filename;
    UniqueThreeVectors<double> Kpoints;
    bool PrintPosition=true, PrintFrequency=true, PrintIntensity=true;
};

#endif /* defined(__SpinWaveDispersion__) */
