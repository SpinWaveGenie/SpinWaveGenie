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
#include "EnergyResolutionFunction.h"
#include "OneDimensionalFactory.h"
#include "OneDimensionalShapes.h"
#include "IntegrateAxes.h"
#include "PointsAlongLine.h"
#include "TwoDimensionCut.h"
#include "Containers/HKLDirections.h"
#include "External/ezRateProgressBar.hpp"


using namespace std;
using namespace tbb;

tbb::atomic<int> counter = 0;
Eigen::MatrixXd mat;
ez::ezRateProgressBar<int> p(1601);

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
            p.update(counter);
        }
    }
    
    
    ApplyFoo(unique_ptr<SpinWavePlot> inCut, ThreeVectors<double> inPoints, size_t inEnergyPoints) : cut(move(inCut)), points(inPoints), energyPoints(inEnergyPoints)
    {
    }
};


void ParallelApplyFoo(unique_ptr<SpinWavePlot> cut, ThreeVectors<double> points, size_t inEnergyPoints)
{
    parallel_for(blocked_range<size_t>(0,points.size()), ApplyFoo(move(cut),points,inEnergyPoints));
}


int main()
{

    tbb::task_scheduler_init init(8);
    double SA = 2.0;
    double theta = M_PI/2.0;
    
    Cell cell;
    cell.setBasisVectors(4.2827,4.2827,6.1103,90.0,90.0,120.0);
    
    Sublattice Spin0;
    string name = "Spin0";
    Spin0.setName(name);
    Spin0.setType("MN2");
    Spin0.setMoment(SA,theta,0.0);
    cell.addSublattice(Spin0);
    cell.addAtom(name,0.0,0.0,0.0);
    cell.addAtom(name,0.0,0.0,0.5);
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactionFactory;
    
    builder.addInteraction(interactionFactory.getExchange("J1",-5.687500,name,name,3.0,3.1));
    builder.addInteraction(interactionFactory.getExchange("J2",4.860082,name,name,4.2,4.3));
    builder.addInteraction(interactionFactory.getExchange("J3",1.781250,name,name,5.2,5.3));
    builder.addInteraction(interactionFactory.getExchange("J4",1.578125,name,name,6.0,6.2));
    builder.addInteraction(interactionFactory.getExchange("J5",-3.188207,name,name,7.40,7.43));
    builder.addInteraction(interactionFactory.getExchange("J6",0.8203125,name,name,7.44,7.49));
    
    SpinWave SW = builder.Create_Element();
    
    OneDimensionalFactory factory;
    auto gauss = factory.getLorentzian(5.0,0.000001);
    
    size_t numberpoints = 1601;
    
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(gauss), SW, 0.0, 120.0, numberpoints));
    
    HKLDirections directions;
    directions.addDirection(0, 0.05);
    directions.addDirection(0.5,-1.0,0.0,0.1);
    directions.addDirection(2, 0.05);
    
    unique_ptr<SpinWavePlot> cut(new IntegrateAxes(move(res),directions,0.0000001));
    
    mat.resize(numberpoints,numberpoints);
    
    PointsAlongLine calculatePoints;
    calculatePoints.setFirstPoint(-2.0, 0.0,-7.0);
    calculatePoints.setFinalPoint(-2.0, 0.0, 3.0);
    calculatePoints.setNumberPoints(numberpoints);
    ThreeVectors<double> points = calculatePoints.getPoints();
    
   
    //TwoDimensionCut twodimcut;
    //twodimcut.setPlotObject(move(cut2));
    //twodimcut.setPoints(points);
    p.units = "Q-points";
    p.start();
    ParallelApplyFoo(move(cut),points,numberpoints);
    //mat = twodimcut.getMatrix();
    std::ofstream file("MnBi_m20L.txt");
    file << mat << endl;

    return 0;
}
