#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include "Genie/SpinWaveBuilder.h"
#include "Genie/Neighbors.h"
#include "Containers/Cell.h"
#include "Interactions/InteractionFactory.h"
#include "External/ezRateProgressBar.hpp"

using namespace std;
using namespace Eigen;
using namespace SpinWaveGenie;

double SA = 1.3;
double SB = 0.33;
Cell cell;
SpinWaveBuilder builder;


void initialize(const std::vector<double> &x, const std::vector<double> &parameters)
{

    
    /*cout << "x[i] = ";
    for (size_t i = 0; i<2;i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    */

    double theta0 = x[0];
    double theta1 = x[1];

    cell.setBasisVectors(8.5,8.5,8.5,90.0,90.0,90.0);
    
    Sublattice Mn0;
    string name = "Mn0";
    Mn0.setName(name);
    Mn0.setType("MN2");
    Mn0.setMoment(SA,0.0,0.0);
    cell.addSublattice(Mn0);
    cell.addAtom(name,0.0,0.0,0.0);
    cell.addAtom(name,0.0,0.5,0.5);
    cell.addAtom(name,0.5,0.0,0.5);
    cell.addAtom(name,0.5,0.5,0.0);
    
    Sublattice Mn1;
    name = "Mn1";
    Mn1.setName(name);
    Mn1.setType("MN2");
    Mn1.setMoment(SA,0.0,0.0);
    cell.addSublattice(Mn1);
    cell.addAtom(name,0.75,0.25,0.75);
    cell.addAtom(name,0.75,0.75,0.25);
    cell.addAtom(name,0.25,0.25,0.25);
    cell.addAtom(name,0.25,0.75,0.75);
    
    Sublattice V0;
    name = "V0";
    V0.setName(name);
    V0.setType("V3");
    if (theta1 < M_PI)
    {
        V0.setMoment(SB,theta1,3.0*M_PI/4.0);
    }
    else
    {
        V0.setMoment(SB,2.0*M_PI-theta1,7.0*M_PI/4.0);
    }
    cell.addSublattice(V0);
    cell.addAtom(name,0.875,0.125,0.375);
    cell.addAtom(name,0.875,0.625,0.875);
    cell.addAtom(name,0.375,0.125,0.875);
    cell.addAtom(name,0.375,0.625,0.375);
    
    Sublattice V1;
    name = "V1";
    V1.setName(name);
    V1.setType("V3");
    if (theta1 < M_PI)
    {
        V1.setMoment(SB,theta1,7.0*M_PI/4.0);
    }
    else
    {
        V1.setMoment(SB,2.0*M_PI-theta1,3.0*M_PI/4.0);
    }
    cell.addSublattice(V1);
    cell.addAtom(name,0.125,0.375,0.875);
    cell.addAtom(name,0.125,0.875,0.375);
    cell.addAtom(name,0.625,0.375,0.375);
    cell.addAtom(name,0.625,0.875,0.875);
    
    Sublattice V2;
    name = "V2";
    V2.setName(name);
    V2.setType("V3");
    if (theta0 < M_PI)
    {
        V2.setMoment(SB,theta0,M_PI/4.0);
    }
    else
    {
        V2.setMoment(SB,2.0*M_PI-theta0,5.0*M_PI/4.0);
    }
    cell.addSublattice(V2);
    cell.addAtom(name,0.375,0.875,0.125);
    cell.addAtom(name,0.375,0.375,0.625);
    cell.addAtom(name,0.875,0.875,0.625);
    cell.addAtom(name,0.875,0.375,0.125);
    
    Sublattice V3;
    name = "V3";
    V3.setName(name);
    V3.setType("V3");
    if (theta0 < M_PI)
    {
        V3.setMoment(SB,theta0,5.0*M_PI/4.0);
    }
    else
    {
        V3.setMoment(SB,2.0*M_PI-theta0,M_PI/4.0);
    }
    cell.addSublattice(V3);
    cell.addAtom(name,0.625,0.625,0.625);
    cell.addAtom(name,0.625,0.125,0.125);
    cell.addAtom(name,0.125,0.625,0.125);
    cell.addAtom(name,0.125,0.125,0.625);
    
    builder.updateCell(cell);
    InteractionFactory interactions;

    builder.addInteraction(interactions.getExchange("Jbb",parameters[1],"V0","V1",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbb",parameters[1],"V2","V3",2.975,3.06));

    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V0","V2",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V0","V3",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V1","V2",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V1","V3",2.975,3.06));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[3],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[3],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[3],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[3],direction,"V3"));

    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V0",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V1",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V2",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V3",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V0",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V1",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V2",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V3",3.48,3.57));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"Mn0"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"Mn1"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"V0"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"V1"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"V2"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[4],zhat,"V3"));
}

double myfunc(const std::vector<double> &x, const std::vector<double> &parameters)
{
    if (x[0] < M_PI)
    {
        cell.getSublattice("V2").setMoment(SB,x[0],M_PI/4.0);
        cell.getSublattice("V3").setMoment(SB,x[0],5.0*M_PI/4.0);
        
    }
    else
    {
        cell.getSublattice("V2").setMoment(SB,2.0*M_PI-x[0],5.0*M_PI/4.0);
        cell.getSublattice("V3").setMoment(SB,2.0*M_PI-x[0],M_PI/4.0);
    }
    
    if (x[1] < M_PI)
    {
        cell.getSublattice("V0").setMoment(SB,x[1],3.0*M_PI/4.0);
        cell.getSublattice("V1").setMoment(SB,x[1],7.0*M_PI/4.0);

    }
    else
    {
        cell.getSublattice("V0").setMoment(SB,2.0*M_PI-x[1],7.0*M_PI/4.0);
        cell.getSublattice("V1").setMoment(SB,2.0*M_PI-x[1],3.0*M_PI/4.0);
    }
    
    builder.updateCell(cell);
    builder.updateInteraction("H", parameters[4]);
    
    return builder.getEnergy();
}

int main()
{
    int samples(801);

    std::vector<double> angles(2.0*M_PI,2);
    vector<double> parameters = {-2.5,-12.3,-12.3,0.0,0.0};
    initialize(angles,parameters);
    Eigen::MatrixXd results(samples,samples);
    for (size_t frame = 0;frame<501;frame++)
    {
        ez::ezRateProgressBar<int> p(samples*samples);
        p.units = "angles";
        p.start();
        for(int i = 0;i<samples;i++)
        {
            for(int j = 0;j<samples;j++)
            {
                angles[0] = 2.0*M_PI*(double)i/(double)(samples-1);
                angles[1] = 2.0*M_PI*(double)j/(double)(samples-1);
                parameters[4] = frame/10.0;
                results(i,j) = myfunc(angles,parameters);
                p.update(i*samples+j);
            }
        }
        p.update(samples*samples);
        std::ofstream file("BField_"+boost::lexical_cast<string>(frame)+".txt");
        file << results << endl;
    }
    return 0;
}
