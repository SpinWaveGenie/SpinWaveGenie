#ifndef __SpinWaveDispersion__
#define __SpinWaveDispersion__

#include <iostream>
#include "SpinWaveGenie/Containers/ThreeVectors.h"
#include "SpinWaveGenie/Genie/SpinWave.h"

class SpinWaveDispersion
{
public:
    SpinWaveDispersion();
    enum class Options {PrintPosition,PrintFrequency,PrintIntensity};
    void setOptions(Options PrintOptions, bool Value);
    void setFilename(std::string name);
    void setPoints(SpinWaveGenie::ThreeVectors<double> pos);
    void setGenie(SpinWaveGenie::SpinWave SW);
    void save();
private:
    SpinWaveGenie::SpinWave Genie;
    std::string Filename;
    SpinWaveGenie::UniqueThreeVectors<double> Kpoints;
    bool PrintPosition, PrintFrequency, PrintIntensity;
};

#endif /* defined(__SpinWaveDispersion__) */
