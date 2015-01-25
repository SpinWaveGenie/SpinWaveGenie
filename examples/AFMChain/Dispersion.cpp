#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Genie/Genie.h"
#include "SpinWaveGenie/Interactions/Interactions.h"
#include "SpinWaveGenie/Plot/Plot.h"

using namespace std;
using namespace SpinWaveGenie;

int main()
{
    Cell cell;
    cell.setBasisVectors(1.0,10.0,10.0,90.0,90.0,90.0);
    
    Sublattice spin0;
    string name0 = "Spin0";
    spin0.setName(name0);
    spin0.setType("NONE");
    spin0.setMoment(1.0,0.0,0.0);
    cell.addSublattice(spin0);
    cell.addAtom(name0,0.0,0.0,0.0);
    
    Sublattice spin1;
    string name1 = "Spin1";
    spin1.setName(name1);
    spin1.setType("NONE");
    spin1.setMoment(1.0,M_PI,0.0);
    cell.addSublattice(spin1);
    cell.addAtom(name1,0.5,0.0,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    Vector3 xhat(1.0,0.0,0.0);
    builder.addInteraction(interactions.getExchange("J",-1.0,name0,name1,0.4,0.6));
    builder.addInteraction(interactions.getAnisotropy("D",0.1,xhat,name0));
    builder.addInteraction(interactions.getAnisotropy("D",0.1,xhat,name1));

    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(3.0,0.0,0.0);
    Line.setNumberPoints(61);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("AFMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    return 0;
}
