#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <memory>
#include <Eigen/Dense>
#include "Genie/SpinWave.h"
#include "Genie/SpinWaveBuilder.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/InteractionFactory.h"
#include "EnergyResolutionFunction.h"
#include "OneDimensionalFactory.h"
#include "OneDimensionalShapes.h"
#include "IntegrateThetaPhi.h"
#include "PointsAlongLine.h"
#include "TwoDimensionCut.h"


using namespace std;

Eigen::MatrixXd mat;


int main()
{
    double SA = 2.0;
    double beta = 22.5*M_PI/180.0;
    double J1 = -1.0;
    double J2 = 2.0*cos(beta)*J1;
    cout << J1 << " " << J2 << endl;
    //double beta = 0.0;
    //double J1 = -1.25;
    //double J2 = -2.5;
    
    Cell cell;
    cell.setBasisVectors(8.5331,5.7892,18.611,90.0,90.0,90.0);
    
    Sublattice Spin1a;
    string name1a = "Spin1a";
    Spin1a.setName(name1a);
    Spin1a.setType("FE3");
    Spin1a.setMoment(SA,M_PI/2.0,M_PI/2.0);
    cell.addSublattice(Spin1a);
    cell.addAtom(name1a,0.71759,0.75,0.56576);
    
    Sublattice Spin1b;
    string name1b = "Spin1b";
    Spin1b.setName(name1b);
    Spin1b.setType("FE3");
    Spin1b.setMoment(SA,M_PI/2.0,3.0*M_PI/2.0);
    cell.addSublattice(Spin1b);
    cell.addAtom(name1b,0.78241,0.25,0.06576);
    
    Sublattice Spin1c;
    string name1c = "Spin1c";
    Spin1c.setName(name1c);
    Spin1c.setType("FE3");
    Spin1c.setMoment(SA,M_PI/2.0,M_PI/2.0);
    cell.addSublattice(Spin1c);
    cell.addAtom(name1c,0.28241,0.25,0.43424);
    
    Sublattice Spin1d;
    string name1d = "Spin1d";
    Spin1d.setName(name1d);
    Spin1d.setType("FE3");
    Spin1d.setMoment(SA,M_PI/2.0,3.0*M_PI/2.0);
    cell.addSublattice(Spin1d);
    cell.addAtom(name1d,0.21759,0.75,0.93424);
    
    Sublattice Spin2a;
    string name2a = "Spin2a";
    Spin2a.setName(name2a);
    Spin2a.setType("FE3");
    Spin2a.setMoment(SA,M_PI/2.0,3.0*M_PI/2.0-beta);
    cell.addSublattice(Spin2a);
    cell.addAtom(name2a,0.0,0.0,0.5);
    
    Sublattice Spin2b;
    string name2b = "Spin2b";
    Spin2b.setName(name2b);
    Spin2b.setType("FE3");
    Spin2b.setMoment(SA,M_PI/2.0,M_PI/2.0-beta);
    cell.addSublattice(Spin2b);
    cell.addAtom(name2b,0.5,0.0,0.0);
    
    Sublattice Spin2c;
    string name2c = "Spin2c";
    Spin2c.setName(name2c);
    Spin2c.setType("FE3");
    Spin2c.setMoment(SA,M_PI/2.0,3.0*M_PI/2.0+beta);
    cell.addSublattice(Spin2c);
    cell.addAtom(name2c,0.0,0.5,0.5);
    
    Sublattice Spin2d;
    string name2d = "Spin2d";
    Spin2d.setName(name2d);
    Spin2d.setType("FE3");
    Spin2d.setMoment(SA,M_PI/2.0,M_PI/2.0+beta);
    cell.addSublattice(Spin2d);
    cell.addAtom(name2d,0.5,0.5,0.0);
    
    /*Neighbors neighborlist;
    neighborlist.findNeighbors(cell,name1d,name2b,2.8,3.1);
    
   
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << " " << dist << endl;
    }*/
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactionFactory;
    
    builder.addInteraction(interactionFactory.getExchange("J1",J1,"Spin2a","Spin2c",2.8,3.0));
    builder.addInteraction(interactionFactory.getExchange("J1",J1,"Spin2b","Spin2d",2.8,3.0));

    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1a","Spin2a",3.0,3.2));
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1a","Spin2c",3.0,3.2));
    
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1c","Spin2a",3.0,3.2));
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1c","Spin2c",3.0,3.2));
    
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1b","Spin2b",3.0,3.2));
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1b","Spin2d",3.0,3.2));
    
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1d","Spin2b",3.0,3.2));
    builder.addInteraction(interactionFactory.getExchange("J2",J2,"Spin1d","Spin2d",3.0,3.2));
    
    builder.addInteraction(interactionFactory.getAnisotropy("D1", -0.1, Vector3(0.0,1.0,0.0), "Spin1a"));
    builder.addInteraction(interactionFactory.getAnisotropy("D1", -0.1, Vector3(0.0,1.0,0.0), "Spin1b"));
    builder.addInteraction(interactionFactory.getAnisotropy("D1", -0.1, Vector3(0.0,1.0,0.0), "Spin1c"));
    builder.addInteraction(interactionFactory.getAnisotropy("D1", -0.1, Vector3(0.0,1.0,0.0), "Spin1d"));
    
    builder.addInteraction(interactionFactory.getAnisotropy("D2",  0.1, Vector3(0.0,0.0,1.0), "Spin2a"));
    builder.addInteraction(interactionFactory.getAnisotropy("D2",  0.1, Vector3(0.0,0.0,1.0), "Spin2b"));
    builder.addInteraction(interactionFactory.getAnisotropy("D2",  0.1, Vector3(0.0,0.0,1.0), "Spin2c"));
    builder.addInteraction(interactionFactory.getAnisotropy("D2",  0.1, Vector3(0.0,0.0,1.0), "Spin2d"));
    
    SpinWave SW = builder.Create_Element();

    /*
    for (int i=0;i<21;i++)
    {
      cout << i/10.0 << endl;
      SW.createMatrix(0.0,0.0,i/10.0);
      SW.calculate();
    
      vector<point> points = SW.getPoints();
      for(auto pt=points.begin();pt!=points.end();pt++)
      {
        cout << pt->frequency << " " << pt->intensity << endl;
      }
    cout << endl;
    }
    */
    
    OneDimensionalFactory factory;
    auto lorentz = factory.getGaussian(1.0,0.001);
    
    size_t numberpoints = 101;
    
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(lorentz), SW, 0.0, 12.0, numberpoints));
    
    unique_ptr<IntegrateThetaPhi> cut(new IntegrateThetaPhi(move(res),0.01));

    unique_ptr<SpinWavePlot> cut2(move(cut));

    PointsAlongLine calculatePoints;
    calculatePoints.setFirstPoint(0.0, 0.0, 0.0);
    calculatePoints.setFinalPoint(0.0, 0.0, 2.0);
    calculatePoints.setNumberPoints(numberpoints);
    ThreeVectors<double> points = calculatePoints.getPoints();
    
    TwoDimensionCut twodimcut;
    twodimcut.setPlotObject(move(cut2));
    twodimcut.setPoints(points);
    mat = twodimcut.getMatrix();
    std::ofstream file("MnBi_H02.txt");
    file << mat << endl;
    
    return 0;
}
