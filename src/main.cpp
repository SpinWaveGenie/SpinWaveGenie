#include <cmath>
#include <nlopt.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "Genie/SpinWave.h"
#include "Initializer.h"
#include "Cell/Cell.h"
#include "Interactions/InteractionFactory.h"
#include "CalculateAngles.h"


using namespace std;
using namespace Eigen;

double SA = 2.5;
double SB = 0.5;

double getCantingAngle(double JAB, double JBB, double JBBP, double DB)
{
    CalculateAngles angles(SA,SB,JAB,JBB,JBBP,DB);
    
    double minf = 0.0;
    std::vector<double> ubang = {M_PI,2.0*M_PI,M_PI,2.0*M_PI,M_PI,2.0*M_PI,M_PI,2.0*M_PI};
    std::vector<double> lbang = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    std::vector<double> thetaphi = {2.44,3.0*M_PI/4.0,2.44,7.0*M_PI/4.0,2.44,M_PI/4.0,2.44,5.0*M_PI/4.0};
    
    nlopt::opt optangle(nlopt::LN_COBYLA,8);
    optangle.set_upper_bounds(ubang);
    optangle.set_lower_bounds(lbang);
    
    optangle.set_ftol_abs(1.0e-13);
    optangle.set_maxeval(5000);
    
    optangle.set_min_objective(CalculateAngles::calc,&angles);
    
    try
    {
        optangle.optimize(thetaphi,minf);
    }
    catch (...)
    {
        cout << "exception occurred" << endl;
    }
    
    return thetaphi[0];
}

double myfunc(const std::vector<double> &thetavector, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    vector<double> x = {-2.2,-15.2,-15.9,0.2};
    
    double theta = thetavector[0];
    
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
    
    InteractionFactory factory;
    
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn0","V0",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn0","V1",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn0","V2",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn0","V3",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn1","V0",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn1","V1",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn1","V2",3.48,3.57));
    builder.addInteraction(factory.getExchange("Jab",x[0],"Mn1","V3",3.48,3.57));
    
    double jbb = -16.0;
    
    builder.addInteraction(factory.getExchange("Jbb",jbb,"V0","V1",2.975,3.06));
    builder.addInteraction(factory.getExchange("Jbb",jbb,"V2","V3",2.975,3.06));

    double jbbp = -1.0*jbb + x[1];
    builder.addInteraction(factory.getExchange("Jbbp",jbbp,"V0","V2",2.975,3.06));
    builder.addInteraction(factory.getExchange("Jbbp",jbbp,"V0","V3",2.975,3.06));
    builder.addInteraction(factory.getExchange("Jbbp",jbbp,"V1","V2",2.975,3.06));
    builder.addInteraction(factory.getExchange("Jbbp",jbbp,"V1","V3",2.975,3.06));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.addInteraction(factory.getAnisotropy("Db",x[2],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.addInteraction(factory.getAnisotropy("Db",x[2],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.addInteraction(factory.getAnisotropy("Db",x[2],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.addInteraction(factory.getAnisotropy("Db",x[2],direction,"V3"));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.addInteraction(factory.getAnisotropy("Da",x[3],zhat,"Mn0"));
    builder.addInteraction(factory.getAnisotropy("Da",x[3],zhat,"Mn1"));
    
    double energy = builder.getEnergy();
    cout << theta << " " << energy << endl;
    cout << builder.getFirstOrderTerms().norm() << endl;
    return energy;
    
}

int main()
{
    double minf = 0.0;
    std::vector<double> ubang = {M_PI};
    std::vector<double> lbang = {0.0};
    std::vector<double> thetaphi = {2.44};
    
    nlopt::opt optangle(nlopt::LN_COBYLA,1);
    optangle.set_upper_bounds(ubang);
    optangle.set_lower_bounds(lbang);
    
    optangle.set_ftol_abs(1.0e-13);
    optangle.set_maxeval(5000);
    
    optangle.set_min_objective(myfunc,nullptr);
    
    try
    {
        optangle.optimize(thetaphi,minf);
    } catch (...)
    {
        cout << "exception occurred" << endl;
    }
    
    cout << thetaphi[0] << endl;

    cout << getCantingAngle(-2.2,-16.0,0.8,-15.9) << endl;

    return 0;
}
