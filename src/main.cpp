#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <memory>
#include <Eigen/Dense>
#include <unistd.h>
#include "tbb/tbb.h"
#include "Genie/SpinWave.h"
#include "Genie/SpinWaveBuilder.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/InteractionFactory.h"
#include "SpinWavePlot/TwoDimensionalGaussian.h"
#include "SpinWavePlot/EnergyResolutionFunction.h"
#include "SpinWavePlot/OneDimensionalFactory.h"
#include "SpinWavePlot/OneDimensionalShapes.h"
#include "SpinWavePlot/IntegrateAxes.h"
#include "Containers/PointsAlongLine.h"
#include "SpinWavePlot/TwoDimensionCut.h"
#include "Containers/HKLDirections.h"
#include "External/ezRateProgressBar.hpp"


using namespace std;
using namespace tbb;
using namespace SpinWaveGenie;


size_t numberpoints = 801;
tbb::atomic<int> counter = 0;
Eigen::MatrixXd mat;

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
                mat(j,i) = data[j];
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

    tbb::task_scheduler_init init(24);
    
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
    info.tol = 1.0e-9;
    
    //size_t EnergyPoints = 17;
    
    Energies energies(0.0,80.0,numberpoints);
    
    //OneDimensionalFactory factory;
    //auto gauss = factory.getLorentzian(5.0,0.000001);
    //unique_ptr<SpinWavePlot> gaussian(new EnergyResolutionFunction(move(gauss),YFeO3,energies));
        
    unique_ptr<SpinWavePlot> gaussian(new TwoDimensionResolutionFunction(info,YFeO3,energies));
    
    HKLDirections directions;
    directions.addDirection(0, 0.2);
    directions.addDirection(1,0.05);
    directions.addDirection(2, 0.2);
    
    unique_ptr<SpinWavePlot> cut(new IntegrateAxes(move(gaussian),directions,1.0e-8));
    
    PointsAlongLine line;
    line.setFirstPoint(2.0,-1.5,-3.0);
    line.setFinalPoint(2.0,1.5,-3.0);
    line.setNumberPoints(numberpoints);
    ThreeVectors<double> points = line.getPoints();

    //TwoDimensionCut twodimcut;
    //twodimcut.setPlotObject(move(cut2));
    //twodimcut.setPoints(points);

    mat.resize(numberpoints,numberpoints);

    tbb::tbb_thread myThread(MyThread);
    ParallelApplyFoo(move(cut),points,numberpoints);
    myThread.join();
    //mat = twodimcut.getMatrix();
    std::ofstream file("MnBi_m20L.txt");
    file << mat << endl;

    return 0;
}
