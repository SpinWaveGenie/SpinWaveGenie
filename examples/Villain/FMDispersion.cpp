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
    cell.setBasisVectors(1.0,2.0,10.0,90.0,90.0,90.0);
    
    Sublattice a1;
    string name1 = "a1";
    a1.setName(name1);
    a1.setType("NONE");
    a1.setMoment(1.0,0.0,0.0);
    cell.addSublattice(a1);
    cell.addAtom(name1,0.0,0.0,0.0);
    
    Sublattice b1;
    string name3 = "b1";
    b1.setName(name3);
    b1.setType("NONE");
    b1.setMoment(1.0,0.0,0.0);
    cell.addSublattice(b1);
    cell.addAtom(name3,0.0,0.5,0.0);
    
    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    double gamma,eta,J,B;
    
    gamma = 3.0;
    eta = -5.0;
    J = 1.0;
    B = 1.0;
    
    builder.addInteraction(interactions.getExchange("J",J,name1,name1,0.8,1.2));
    builder.addInteraction(interactions.getExchange("metaJ",-eta*J,name3,name3,0.8,1.2));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name1,name3,0.8,1.2));
    
    //Vector3 zhat(0.0,0.0,1.0);
    //builder.addInteraction(interactions.getMagneticField("B",B/2.0,zhat,name1));
    //builder.addInteraction(interactions.getMagneticField("B",B/2.0,zhat,name3));


    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(1.0,0.0,0.0);
    Line.setFinalPoint(2.0,0.0,0.0);
    Line.setNumberPoints(11);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("AFMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    
    
    double kx(0.0),ky(0.0);
    double R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    double term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    
    kx = M_PI;
    R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    
    kx = 2.0*M_PI;
    R1K = sqrt(pow(eta+1,2)*pow(1-cos(kx),2)+4.0*pow(gamma*cos(ky),2));
    term = B + 2*gamma + (eta-1.0)*(cos(kx)-1.0);
    cout << term+R1K << " " << term-R1K << endl;
    cout << (R1K-2*gamma*cos(ky))/R1K << " " << (R1K+2*gamma*cos(ky))/R1K << endl;

    return 0;
}
