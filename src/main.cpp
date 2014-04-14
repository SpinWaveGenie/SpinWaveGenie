#include <cmath>
#include <nlopt.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Genie/SpinWave.h"
#include "Initializer.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/ExchangeInteraction.h"
#include "Interactions/AnisotropyInteraction.h"
#include "TwoDimensionCut.h"
#include "OneDimensionalFactory.h"
#include "OneDimensionalShapes.h"
#include "PointsAlongLine.h"
#include "CalculateAngles.h"

using namespace std;
using namespace Eigen;


double myfunc(const std::vector<double> &x)
{

    cout << "x[i] = ";
    for (size_t i = 0; i< x.size();i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    
    double SA = 2.5;
    double SB = 2.5/4.01;
    
    double JAB = x[0];
    double JBB = -7.73362;
    double JBBP = -1.8592;
    double DB = x[1];
    double DA = 0.167;

    CalculateAngles angles(SA,SB,JAB,JBB,JBBP,DB);
    
    double minf = 0.0;
    std::vector<double> ub = {M_PI,2.0*M_PI,M_PI,2.0*M_PI,M_PI,2.0*M_PI,M_PI,2.0*M_PI};
    std::vector<double> lb = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    std::vector<double> thetaphi = {2.44,3.0*M_PI/4.0,2.44,7.0*M_PI/4.0,2.44,M_PI/4.0,2.44,5.0*M_PI/4.0};

    nlopt::opt opt(nlopt::LN_COBYLA,8);
    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);
    
    opt.set_ftol_abs(1.0e-12);
    opt.set_maxeval(5000);
    
    opt.set_min_objective(CalculateAngles::calc,&angles);
    opt.optimize(thetaphi,minf);
    
    
    double theta = thetaphi[0];
    
    cout << "theta = " << 180/M_PI*theta << endl;

    Cell cell;
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
    V0.setMoment(SB,theta,3.0*M_PI/4.0);
    cell.addSublattice(V0);
    cell.addAtom(name,0.875,0.125,0.375);
    cell.addAtom(name,0.875,0.625,0.875);
    cell.addAtom(name,0.375,0.125,0.875);
    cell.addAtom(name,0.375,0.625,0.375);
    
    Sublattice V1;
    name = "V1";
    V1.setName(name);
    V1.setType("V3");
    V1.setMoment(SB,theta,7.0*M_PI/4.0);
    cell.addSublattice(V1);
    cell.addAtom(name,0.125,0.375,0.875);
    cell.addAtom(name,0.125,0.875,0.375);
    cell.addAtom(name,0.625,0.375,0.375);
    cell.addAtom(name,0.625,0.875,0.875);
    
    Sublattice V2;
    name = "V2";
    V2.setName(name);
    V2.setType("V3");
    V2.setMoment(SB,theta,M_PI/4.0);
    cell.addSublattice(V2);
    cell.addAtom(name,0.375,0.875,0.125);
    cell.addAtom(name,0.375,0.375,0.625);
    cell.addAtom(name,0.875,0.875,0.625);
    cell.addAtom(name,0.875,0.375,0.125);
    
    Sublattice V3;
    name = "V3";
    V3.setName(name);
    V3.setType("V3");
    V3.setMoment(SB,theta,5.0*M_PI/4.0);
    cell.addSublattice(V3);
    cell.addAtom(name,0.625,0.625,0.625);
    cell.addAtom(name,0.625,0.125,0.125);
    cell.addAtom(name,0.125,0.625,0.125);
    cell.addAtom(name,0.125,0.125,0.625);
    
    SW_Builder builder(cell);
    
    builder.Add_Interaction(new ExchangeInteraction("Jbb",JBB,"V0","V1",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbb",JBB,"V2","V3",2.975,3.06));

    builder.Add_Interaction(new ExchangeInteraction("Jbbp",JBBP,"V0","V2",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",JBBP,"V0","V3",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",JBBP,"V1","V2",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",JBBP,"V1","V3",2.975,3.06));
       
    Vector3 direction(-1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",DB,direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",DB,direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",DB,direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",DB,direction,"V3"));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Daz",DA,zhat,"Mn0"));
    builder.Add_Interaction(new AnisotropyInteraction("Daz",DA,zhat,"Mn1"));
    
    //double JBB = x[0];
    //double JBBP = x[1];
    //double DB = x[2];
    //double JAB = SB*((6.0*JBB+6.0*JBBP+DB)*cos(theta)*sin(theta) -sqrt(2.0)*DB*(2.0*pow(cos(theta),2)-1.0))/(-9.0*SA*sin(theta));
    
    cout << "JAB= " << JAB << endl;
    
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn0","V0",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn0","V1",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn0","V2",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn0","V3",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn1","V0",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn1","V1",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn1","V2",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",JAB,"Mn1","V3",3.48,3.57));
    
    SpinWave SW = builder.Create_Element();
    SW.createMatrix(0.0,0.0,0.0);
    SW.Calc();
    vector<point> points = SW.getPoints();
    return points[0].frequency;
}

int main()
{
    cout << "#MinValue JAB JBB JBBP DAZ DBZ DB111 " << endl;

    std::vector<double> ub(2),lb(2);
    std::vector<size_t> ns(2);
    ub[0] =  -1.5;
    ub[1] =  -8.0;

    lb[0] = -2.5;
    lb[1] = -10.0;
    
    ns[0] =  100;
    ns[1] =  100;
    
    
    MatrixXd plot;
    plot.resize(ns[0]+1.0,ns[1]+1.0);
    
    for (int i=0;i<=ns[0];i++)
    {
        double jab = lb[0]*(1.0-(double)i/(double)ns[0]) + ub[0]*(double)i/(double)ns[0];
        for (int j=0;j<=ns[1];j++)
        {
            double db = lb[1]*(1.0-(double)j/(double)ns[1]) + ub[1]*(double)j/(double)ns[1];
            vector<double> xtemp = {jab,db};
            plot(i,j) = myfunc(xtemp);
        }
    }
    
    std::ofstream file("jab_db.txt");
    if (file.is_open())
    {
        file << plot;
    }
    file << endl;
    file.close();
    
    cout << plot << endl << endl;
    cout << plot(0,0) << " " << plot(0,200) << " " << plot(200,0) << " " << plot(200,200) << endl;

    return 0;
}
