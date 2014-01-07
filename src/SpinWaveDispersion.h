#ifndef __SpinWaveDispersion__
#define __SpinWaveDispersion__

#include <iostream>
#include "Containers/Positions.h"
#include "Genie/SpinWave.h"

class SpinWaveDispersion
{
public:
    enum class Options {PrintPosition,PrintFrequency,PrintIntensity};
    void setOptions(Options PrintOptions, bool Value);
    void setFilename(std::string name);
    void setPoints(Positions pos);
    void setGenie(SpinWave SW);
    void save();
private:
    SpinWave Genie;
    std::string Filename;
    Positions Kpoints;
    bool PrintPosition=true, PrintFrequency=true, PrintIntensity=true;
};

#endif /* defined(__SpinWaveDispersion__) */
