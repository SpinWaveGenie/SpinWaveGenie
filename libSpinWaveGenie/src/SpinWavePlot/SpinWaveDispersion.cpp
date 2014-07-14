#include <fstream>
#include "SpinWaveGenie/SpinWavePlot/SpinWaveDispersion.h"
#include "SpinWaveGenie/Containers/Results.h"

using std::string; using std::vector;
using std::cout; using std::endl;
using namespace SpinWaveGenie;


SpinWaveDispersion::SpinWaveDispersion()
{
    PrintPosition = true;
    PrintFrequency = true;
    PrintIntensity = true;
}

void SpinWaveDispersion::setOptions(Options PrintOptions, bool Value)
{
    switch(PrintOptions)
    {
        case (Options::PrintPosition):
            PrintPosition = Value;
        case (Options::PrintFrequency):
            PrintFrequency = Value;
        case (Options::PrintIntensity):
            PrintIntensity = Value;
    }
}

void SpinWaveDispersion::setFilename(string name)
{
    Filename = name;
}

void SpinWaveDispersion::setGenie(SpinWave SW)
{
    Genie = SW;
}

void SpinWaveDispersion::setPoints(ThreeVectors<double> pts)
{
    for(auto it = pts.begin(); it!= pts.end(); it++)
    {
        Kpoints.insert(it->get<0>(),it->get<1>(),it->get<2>());
    }
}

void SpinWaveDispersion::save()
{
    
    std::ofstream file(Filename);
    for(auto it = Kpoints.begin(); it != Kpoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();
        cout << x << " " << y << " " << z << endl;

        if (PrintPosition)
        {
            file << x << " " << y << " " << z << " ";// << endl;
        }
        
        Genie.createMatrix(x,y,z);
        Genie.calculate();
        Results pts = Genie.getPoints();
        
        if(PrintFrequency)
        {
            for(Results::Iterator it = pts.begin();it!=pts.end();it++)
            {
                file << it->frequency << "  ";
            }
        }
        
        if(PrintIntensity)
        {
            for(Results::Iterator it = pts.begin();it!=pts.end();it++)
            {
                file << it->intensity << " " ;
            }
        }
        
        file << endl;
    }
}


