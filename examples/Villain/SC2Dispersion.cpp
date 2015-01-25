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
    
    double gamma,eta,J,B;
    
    gamma = -2.0;
    eta = -2.0;
    J = 1.0;
    B = 8.0;
    
    double theta = M_PI_2;
    theta = acos(-B/(4.0*gamma));
    cout << "Theta = " << theta << endl;
    
    Cell cell;
    cell.setBasisVectors(1.0,2.0,10.0,90.0,90.0,90.0);
    
    Sublattice a1;
    string name1 = "a1";
    a1.setName(name1);
    a1.setType("NONE");
    a1.setMoment(1.0,theta,0.0);
    cell.addSublattice(a1);
    cell.addAtom(name1,0.0,0.0,0.0);
    
    Sublattice b1;
    string name2 = "b1";
    b1.setName(name2);
    b1.setType("NONE");
    b1.setMoment(1.0,theta,M_PI);
    cell.addSublattice(b1);
    cell.addAtom(name2,0.0,0.5,0.0);

    SpinWaveBuilder builder(cell);
    
    InteractionFactory interactions;
    
    builder.addInteraction(interactions.getExchange("J",J,name1,name1,0.9,1.1));
    builder.addInteraction(interactions.getExchange("metaJ",-1.0*eta*J,name2,name2,0.9,1.1));
    builder.addInteraction(interactions.getExchange("gammaJ",gamma*J,name1,name2,0.9,1.1));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name1));
    builder.addInteraction(interactions.getMagneticField("B",B,zhat,name2));

    SpinWave test = builder.createElement();
    
    PointsAlongLine Line;
    Line.setFirstPoint(0.0,1.0,0.0);
    Line.setFinalPoint(0.0,2.0,0.0);
    Line.setNumberPoints(11);
    ThreeVectors<double> kPoints = Line.getPoints();
    
    SpinWaveDispersion dispersion;
    dispersion.setFilename("AFMChain.txt");
    dispersion.setGenie(test);
    dispersion.setPoints(kPoints);
    
    dispersion.save();
    
    /*double kx(0.0),ky(0.0);
    double R2Kp = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0+cos(ky));
    double R2Km = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0-cos(ky));
    double term = (eta+1.0)*(cos(kx)-1.0);
    cout << sqrt(R2Kp*R2Km) +term << " " << sqrt(R2Kp*R2Km) +term << endl;
    cout << sqrt(R2Kp/R2Km) << " " << sqrt(R2Kp/R2Km) << endl;
    
    kx = M_PI;
    R2Kp = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0+cos(ky));
    R2Km = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0-cos(ky));
    term = (eta+1.0)*(cos(kx)-1.0);
    cout << sqrt(R2Kp*R2Km) +term << " " << sqrt(R2Kp*R2Km) -term << endl;
    cout << sqrt(R2Kp/R2Km) << " " << sqrt(R2Kp/R2Km) << endl;

    kx = 2.0*M_PI;
    R2Kp = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0+cos(ky));
    R2Km = (eta-1.0)*(cos(kx)-1.0)+2.0*gamma*(-1.0-cos(ky));
    term = (eta+1.0)*(cos(kx)-1.0);
    cout << sqrt(R2Kp*R2Km) +term << " " << sqrt(R2Kp*R2Km) +term << endl;
    cout << sqrt(R2Kp/R2Km) << " " << sqrt(R2Kp/R2Km) << endl;
*/
    return 0;
}
