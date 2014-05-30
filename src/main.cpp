#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
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
#include <unistd.h>
#include "tbb/tbb.h"
#include "External/ezRateProgressBar.hpp"
#include "nexus/NeXusFile.hpp"
#include "Containers/Energies.h"

using namespace std;
using namespace tbb;

size_t numberpoints = 1601;
tbb::atomic<int> counter;
Eigen::MatrixXd mat(numberpoints,numberpoints);

class ApplyFoo
{
    unique_ptr<SpinWavePlot> cut;
    ThreeVectors<double> points;
    size_t energyPoints;
public:
    ApplyFoo(const ApplyFoo& other)
    {
        cut = move(other.cut->clone());
        points = other.points;
        energyPoints = other.energyPoints;
    }
    void operator()( const blocked_range<size_t>& r ) const
    {
        unique_ptr<SpinWavePlot> cut2(move(cut->clone()));
        auto it = points.cbegin()+ r.begin();
        for(size_t i = r.begin(); i!=r.end(); ++i)
        {
            vector<double> data = cut2->getCut(it->get<0>(),it->get<1>(),it->get<2>());
            for (auto j=0;j!=energyPoints;j++)
            {
                mat(i,j) = data[j];
            }
            it++;
            counter++;
        }
    }
    
    
    ApplyFoo(unique_ptr<SpinWavePlot> inCut, ThreeVectors<double> inPoints, size_t inEnergyPoints) : cut(move(inCut)), points(inPoints), energyPoints(inEnergyPoints)
    {
    }
};

void MyThread()
{
    ez::ezRateProgressBar<int> p(numberpoints);
    p.units = "Q-points";
    p.start();
    while( counter < numberpoints)
    {
        p.update(counter);
        sleep(1);
    }
    p.update(numberpoints);
}

void ParallelApplyFoo(unique_ptr<SpinWavePlot> cut, ThreeVectors<double> points, size_t inEnergyPoints)
{
    parallel_for(blocked_range<size_t>(0,points.size()), ApplyFoo(move(cut),points,inEnergyPoints));
}


int main()
{
    tbb::task_scheduler_init init(8);
    double SA = 1.5;
    
    Cell cell;
    cell.setBasisVectors(9.91510,8.84410,5.45950,90.0,107.5780,90.0);
    
    Sublattice Spin0;
    string name0 = "Spin0";
    Spin0.setName(name0);
    Spin0.setType("CR3");
    Spin0.setMoment(SA,17.568*M_PI/180.0,M_PI);
    cell.addSublattice(Spin0);
    cell.addAtom(name0,0.0,0.91226,0.25);
    cell.addAtom(name0,0.5,0.41226,0.25);
    
    Sublattice Spin1;
    string name1 = "Spin1";
    Spin1.setName(name1);
    Spin1.setType("CR3");
    Spin1.setMoment(SA,17.568*M_PI/180.0,M_PI);
    cell.addSublattice(Spin1);
    cell.addAtom(name1,0.0,0.08774,0.75);
    cell.addAtom(name1,0.5,0.58774,0.75);
    
    /*Neighbors neighborlist;
    neighborlist.findNeighbors(cell, name0, name1, 5.6, 5.7);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << " " << dist << endl;
    }
    */
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactionFactory;
    
    builder.addInteraction(interactionFactory.getExchange("J",0.448,name0,name1,3.1,3.2));
    builder.addInteraction(interactionFactory.getExchange("J2",0.05,name0,name1,5.6,5.7));
    builder.addInteraction(interactionFactory.getExchange("J1",-0.01,name0,name0,6.6,6.7));
    builder.addInteraction(interactionFactory.getExchange("J1",-0.01,name1,name1,6.6,6.7));
    
    builder.addInteraction(interactionFactory.getAnisotropy("D", -0.1, Vector3(0.0,sin(17.568*M_PI/180.0),-cos(17.568*M_PI/180.0)),"Spin0"));
    builder.addInteraction(interactionFactory.getAnisotropy("D", -0.1, Vector3(0.0,sin(17.568*M_PI/180.0),-cos(17.568*M_PI/180.0)),"Spin1"));

    SpinWave SW = builder.Create_Element();

    /*SW.createMatrix(0.0, 0.0, 0.9);
    SW.Calc();
    
    vector<point> points = SW.getPoints();
    for(auto pt=points.begin();pt!=points.end();pt++)
    {
        cout << pt->frequency << " " << pt->intensity << endl;
    }*/
    
    //TwoDimensionCut twodimcut;
    //twodimcut.setPlotObject(move(cut));
    //twodimcut.setPoints(points);
    //mat = twodimcut.getMatrix();
    
    
    OneDimensionalFactory factory;
    auto lorentz = factory.getGaussian(0.3,0.000001);
    
    Energies energies(0.0, 6.0, numberpoints);
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(lorentz), SW, energies));
    
    unique_ptr<SpinWavePlot> cut(new IntegrateThetaPhi(move(res),0.001));

    PointsAlongLine calculatePoints;
    calculatePoints.setFirstPoint( 0.0, 0.0, 0.0);
    calculatePoints.setFinalPoint( 0.0, 0.0, 2.0);
    calculatePoints.setNumberPoints(numberpoints);
    ThreeVectors<double> points = calculatePoints.getPoints();
    
    tbb::tbb_thread myThread(MyThread);
    ParallelApplyFoo(move(cut),points,numberpoints);
    myThread.join();
    
    //std::ofstream file("MnBi_m20L.txt");
    //file << mat << endl;
    

    NeXus::File nf("MnBi_m20L.hdf5",NXACC_CREATE5);
    nf.makeGroup("entry","NXentry",true);
    nf.makeGroup("data","NXdata",true);
    
    std::vector<int> array_dims;
    array_dims.push_back(numberpoints);
    array_dims.push_back(numberpoints);
    nf.makeData("data", NeXus::FLOAT64, array_dims, true);
    nf.putData(mat.data());
    nf.putAttr("signal", 1);
    nf.putAttr("axes","wavevector:energy");
    nf.putAttr("units","arbitrary units");
    
    nf.closeData();
    
    vector<double> hpts,kpts,lpts;
    for (auto it=points.begin();it!=points.end();++it)
    {
        hpts.push_back(it->get<0>());
        kpts.push_back(it->get<1>());
        lpts.push_back(it->get<2>());
    }
    
    nf.makeData("origin", NeXus::FLOAT64, 4 , true);
    vector<double> origin = {1.0,0.0,-2.0,0.0};
    nf.putData(&(origin[0]));
    nf.closeData();
   
    array_dims.clear();
    array_dims.push_back(2);
    array_dims.push_back(4);
    nf.makeData("axes_index", NeXus::INT32, array_dims , true);
    array<int,4*2> axes_index = {{0,0,1,0, 0,0,0,1}};
    nf.putData(axes_index.data());
    nf.closeData();

    nf.makeData("L", NeXus::FLOAT64, (int)numberpoints , true);
    nf.putData(lpts.data());
    nf.putAttr("axis",1);
    nf.putAttr("primary",1);
    nf.putAttr("units","NX_UNITLESS");
    nf.closeData();
    
    nf.makeData("energies", NeXus::FLOAT64, (int)numberpoints, true);
    nf.putData(energies.data());
    nf.putAttr("axis",2);
    nf.putAttr("units","meV");
    
    nf.closeGroup();
    nf.closeGroup();
    nf.closeData();

    return 0;
}
