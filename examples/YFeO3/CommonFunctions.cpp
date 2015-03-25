#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "CommonFunctions.h"

using namespace SpinWaveGenie;

SpinWaveGenie::SpinWave createModel()
{
    Cell cell;
    cell.setBasisVectors(5.4,5.4,7.63675323681,90.0,90.0,90.0);
    
    double delta = 0.00516;
    double phi = 0.0032;
    
    Sublattice Fe1;
    std::string name1 = "Fe1";
    Fe1.setName(name1);
    Fe1.setType("FE3");
    Fe1.setMoment(2.5,M_PI/2.0-delta,M_PI+phi);
    cell.addSublattice(Fe1);
    cell.addAtom(name1,0.0,0.5,0.0);

    Sublattice Fe2;
    std::string name2 = "Fe2";
    Fe2.setName(name2);
    Fe2.setType("FE3");
    Fe2.setMoment(2.5,M_PI/2.0-delta,phi);
    cell.addSublattice(Fe2);
    cell.addAtom(name2,0.0,0.5,0.5);

    Sublattice Fe3;
    std::string name3 = "Fe3";
    Fe3.setName(name3);
    Fe3.setType("FE3");
    Fe3.setMoment(2.5,M_PI/2.0-delta,M_PI-phi);
    cell.addSublattice(Fe3);
    cell.addAtom(name3,0.5,0.0,0.5);

    Sublattice Fe4;
    std::string name4 = "Fe4";
    Fe4.setName(name4);
    Fe4.setType("FE3");
    Fe4.setMoment(2.5,M_PI/2.0-delta,2.0*M_PI-phi);
    cell.addSublattice(Fe4);
    cell.addAtom(name4,0.5,0.0,0.0);
    
    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
   
    Vector3 xhat(1.0,0.0,0.0); 
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name1));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name2));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name3));
    builder.addInteraction(interactions.getAnisotropy("Ka",-0.0055,xhat,name4));

    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name1));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name2));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name3));
    builder.addInteraction(interactions.getAnisotropy("Kc",-0.00305,zhat,name4));

    builder.addInteraction(interactions.getExchange("J1",-4.77,name1,name2,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name1,name4,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name3,name2,3.8,4.3));
    builder.addInteraction(interactions.getExchange("J1",-4.77,name3,name4,3.8,4.3));

    builder.addInteraction(interactions.getExchange("J2",-0.21,name1,name1,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name2,name2,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name3,name3,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name4,name4,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name1,name3,5.3,5.5));
    builder.addInteraction(interactions.getExchange("J2",-0.21,name2,name4,5.3,5.5));

    Vector3 yhat(0.0,1.0,0.0);
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.074,yhat,name4,name1,3.8,4.3));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D1",-0.074,yhat,name2,name3,3.8,4.3));

    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.028,zhat,name4,name1,3.8,4.3));
    builder.addInteraction(interactions.getDzyaloshinskiiMoriya("D2",-0.028,zhat,name3,name2,3.8,4.3));

    return builder.createElement();
}
