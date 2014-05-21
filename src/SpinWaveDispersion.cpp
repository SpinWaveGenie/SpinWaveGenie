#include "SpinWaveDispersion.h"
#include <fstream>

using std::string; using std::vector;
using std::cout; using std::endl;

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
        vector<point> pts = Genie.getPoints();
        
        if(PrintFrequency)
        {
            for(vector<point>::iterator it = pts.begin();it!=pts.end();it++)
            {
                file << (*it).frequency << "  ";
            }
        }
        
        if(PrintIntensity)
        {
            for(vector<point>::iterator it = pts.begin();it!=pts.end();it++)
            {
                file << (*it).intensity << " " ;
            }
        }
        
        file << endl;
    }
}


