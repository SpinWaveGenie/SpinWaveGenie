//
//  main.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//
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

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",1.0,name0,name0,0.9,1.1));

    SpinWave SW = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,0.0,0.0);
    Line.setFinalPoint(0.0,0.0,3.0*2.0*M_PI);
    Line.setNumberPoints(401);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    Energies energies(0.0, 5.0, 401);
    
    OneDimensionalFactory factory;
    auto gauss = factory.getGaussian(0.25,1.0e-2);
    
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(gauss), SW,energies));
    unique_ptr<SpinWavePlot> cut(new IntegrateThetaPhi(move(res),2.0e-2));

    TwoDimensionalCut twodimcut;
    twodimcut.setFilename("FMPowderAverage");
    twodimcut.setPlotObject(move(cut));
    twodimcut.setPoints(kPoints);
    twodimcut.save();
    return 0;
}
