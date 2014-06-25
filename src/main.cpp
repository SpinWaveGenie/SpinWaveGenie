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
#include <vector>
#include "Cell/Cell.h"
#include "Genie/SpinWaveBuilder.h"
#include "Interactions/InteractionFactory.h"
#include "Containers/Results.h"
#include "SpinWavePlot/TwoDGaussian.h"
#include "SpinWavePlot/TwoDimensionalGaussian.h"
#include "SpinWavePlot/IntegrateAxes.h"
#include "External/ezRateProgressBar.hpp"
#include "Containers/PointsAlongLine.h"
#include "Containers/ThreeVectors.h"
#include "Containers/HKLDirections.h"

using namespace std;
using namespace SpinWaveGenie;

int main(int argc, char * argv[])
{
    double S = 2.5;
    double theta = 89.7043537777*M_PI/180.0;
    double delta = 0.18334649444;
    
    Cell cell;
    cell.setBasisVectors(5.4,5.4,7.63675323681,90.0,90.0,90.0);
    
    Sublattice Fe1;
    Fe1.setName("Fe1");
    Fe1.setType("FE3");
    Fe1.setMoment(S,theta,(180.0 + delta)*M_PI/180.0);
    cell.addSublattice(Fe1);
    cell.addAtom("Fe1",0.0,0.5,0.0);
    
    Sublattice Fe2;
    Fe2.setName("Fe2");
    Fe2.setType("FE3");
    Fe2.setMoment(S,theta,delta*M_PI/180.0);
    cell.addSublattice(Fe2);
    cell.addAtom("Fe2",0.0,0.5,0.5);

    Sublattice Fe3;
    Fe3.setName("Fe3");
    Fe3.setType("FE3");
    Fe3.setMoment(S,theta,(180.0-delta)*M_PI/180.0);
    cell.addSublattice(Fe3);
    cell.addAtom("Fe3",0.5,0.0,0.5);
    
    Sublattice Fe4;
    Fe4.setName("Fe4");
    Fe4.setType("FE3");
    Fe4.setMoment(S,theta,(360.0-delta)*M_PI/180.0);
    cell.addSublattice(Fe4);
    cell.addAtom("Fe4",0.5,0.0,0.0);
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J1",-4.8982,"Fe1","Fe2",3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.8982,"Fe1","Fe4",3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.8982,"Fe3","Fe2",3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.8982,"Fe3","Fe4",3.8,4.3));
    
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe1","Fe1",5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe2","Fe2",5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe3","Fe3",5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe4","Fe4",5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe1","Fe3",5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.252455,"Fe2","Fe4",5.3,5.5));
    
    Vector3 xhat(1.0,0.0,0.0);
    Vector3 yhat(0.0,1.0,0.0);
    Vector3 zhat(0.0,0.0,1.0);
    
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.00527847,xhat,"Fe1"));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.00527847,xhat,"Fe2"));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.00527847,xhat,"Fe3"));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.00527847,xhat,"Fe4"));
    
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00294114,zhat,"Fe1"));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00294114,zhat,"Fe2"));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00294114,zhat,"Fe3"));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00294114,zhat,"Fe4"));
    
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.0758301,yhat,"Fe4","Fe1",3.78,4.32));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.0758301,yhat,"Fe2","Fe3",3.78,4.32));
    
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.0281255,zhat,"Fe4","Fe1",3.78,4.32));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.0281255,zhat,"Fe3","Fe2",3.78,4.32));

    SpinWave YFeO3 = builder.Create_Element();

    TwoDimGaussian info;
    info.a = 1109.0;
    info.b = 0.0;
    info.c = 0.48;
    info.direction = Vector3(0.0,1.0,0.0);
    info.tol = 1.0e-6;
    
    size_t EnergyPoints = 17;
    
    Energies energies(0.0,80.0,EnergyPoints);
    
    unique_ptr<SpinWavePlot> gaussian(new TwoDimensionResolutionFunction(info,YFeO3,energies));
    
    HKLDirections directions;
    directions.addDirection(0, 0.2);
    //directions.addDirection(1,0.05);
    directions.addDirection(2, 0.2);
    
    unique_ptr<SpinWavePlot> cut(new IntegrateAxes(move(gaussian),directions,1.0e-4));
    
    PointsAlongLine line;
    line.setFirstPoint(2.0,-1.5,-3.0);
    line.setFinalPoint(2.0,1.5,-3.0);
    line.setNumberPoints(101);
    ThreeVectors<double> Kpoints = line.getPoints();
    
    Eigen::MatrixXd figure;
    figure.setZero(Kpoints.size(),EnergyPoints);
    ez::ezRateProgressBar<int> p(Kpoints.size());
    p.units = "Q-points";
    p.start();
    for(auto it = Kpoints.begin(); it != Kpoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();
        
        vector<double> val = cut->getCut(x,y,z);
        size_t m = std::distance(Kpoints.begin(),it);
        p.update(m);
        for(int n=0;n<EnergyPoints;n++)
        {
            if (val[n] < 0.0 && std::abs(val[n]) > 1.0e-4)
                cout << "0.0 > val[" << n << "] = " << val[n] << endl;
            figure(m,n) = val[n];
        }
    }
    p.update(Kpoints.size());
 
    
    cout << figure << endl;
    
    /*for (int i=0;i<figure.size();i++)
    {
        cout << i << " " << figure.row(i).sum() << endl;
    }*/
    
    //cout << figure.colwise().sum() << endl;
    //cout << figure.rowwise().sum() << endl;

    return 0;
}
