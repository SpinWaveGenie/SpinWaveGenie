#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <memory>
#include "Genie/SpinWave.h"
#include "Genie/SW_Builder.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/InteractionFactory.h"
#include "EnergyResolutionFunction.h"
#include "OneDimensionalFactory.h"
#include "OneDimensionalShapes.h"
#include "IntegrateAxes.h"
#include <Eigen/Dense>

using namespace std;

int main()
{
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
    
    SW_Builder builder(cell);
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
    
    IntegrateAxes cut(move(res),0.0000001);
    cut.addDirection(0, 0.05);
    cut.addDirection(0.5,-1.0,0.0,0.1);
    cut.addDirection(2, 0.2);

    Eigen::MatrixXd mat;
    mat.resize(numberpoints,numberpoints);
    
    for (int i = 0; i!=numberpoints; i++)
    {
        double x = -4.0+i*7.0/(numberpoints-1);
        cout << x << endl;
        vector<double> data = cut.getCut(x, 0.0 ,2.0);
        for (auto j=0;j!=numberpoints;j++)
        {
            mat(j,i) = data[j];
        }
    }
    
    std::ofstream file("MnBi_H02.txt");
    file << mat << endl;

    return 0;
}
