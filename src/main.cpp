#include <cmath>
#include <nlopt.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Genie/SpinWaveBuilder.h"
#include "Initializer.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/InteractionFactory.h"
//#include "CalculateAngles.h"


using namespace std;
using namespace Eigen;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    double* tmp = static_cast<double*>(my_func_data);
    vector<double> parameters(tmp,tmp+6);
    
    /*cout << "x[i] = ";
    for (size_t i = 0; i<6;i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    */
    
    double SA = 1.3;
    double SB = 0.25;
    double theta0 = x[0];
    double theta1 = x[1];
    
    
    //CalculateAngles angles(SA,SB,parameters[0],parameters[1],parameters[2],parameters[4]);
    
    //vector<double> anglesIn = {theta,3.0*M_PI/4.0,theta,7.0*M_PI/4.0,theta,M_PI/4.0,theta,5.0*M_PI/4.0};
    
    //cout << angles.calculateEnergy(anglesIn) << endl;
    
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
    if (theta1 < M_PI)
    {
        V2.setMoment(SB,theta0,5.0*M_PI/4.0);
    }
    else
    {
        V2.setMoment(SB,2.0*M_PI-theta0,M_PI/4.0);
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
    if (theta1 < M_PI)
    {
        V3.setMoment(SB,theta0,M_PI/4.0);
    }
    else
    {
        V3.setMoment(SB,2.0*M_PI-theta0,5.0*M_PI/4.0);
    }
    cell.addSublattice(V3);
    cell.addAtom(name,0.625,0.625,0.625);
    cell.addAtom(name,0.625,0.125,0.125);
    cell.addAtom(name,0.125,0.625,0.125);
    cell.addAtom(name,0.125,0.125,0.625);
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactions;

    builder.addInteraction(interactions.getExchange("Jbb",parameters[1],"V0","V1",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbb",parameters[1],"V2","V3",2.975,3.06));

    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V0","V2",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V0","V3",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V1","V2",2.975,3.06));
    builder.addInteraction(interactions.getExchange("Jbbp",parameters[2],"V1","V3",2.975,3.06));
    
    Vector3 zhat(0.0,0.0,1.0);

    builder.addInteraction(interactions.getAnisotropy("Daz",parameters[3],zhat,"Mn0"));
    builder.addInteraction(interactions.getAnisotropy("Daz",parameters[3],zhat,"Mn1"));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[4],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[4],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[4],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.addInteraction(interactions.getAnisotropy("Db",parameters[4],direction,"V3"));

    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V0",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V1",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V2",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn0","V3",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V0",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V1",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V2",3.48,3.57));
    builder.addInteraction(interactions.getExchange("Jab",parameters[0],"Mn1","V3",3.48,3.57));
    
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"Mn0"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"Mn1"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"V0"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"V1"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"V2"));
    builder.addInteraction(interactions.getMagneticField("H",parameters[5],zhat,"V3"));
    
    double diff_sq = builder.getEnergy();
    //cout << diff_sq << endl << endl;
    return diff_sq;
}

int main()
{
    /*nlopt::opt opt(nlopt::LN_SBPLX,2);
    std::vector<double> ub(2);
    ub[0] = 2.0*M_PI;
    ub[1] = 2.0*M_PI;
    opt.set_upper_bounds(ub);
    std::vector<double> lb(2);
    lb[0] = 0.0;
    lb[0] = 0.0;
    opt.set_lower_bounds(lb);

    double parameters[6] = {-2.5,-12.0,-12.0,0.0,-1.2,-0.5};
    
    opt.set_ftol_rel(1.0e-8);
    //opt.set_maxeval(1000);

    std::vector<double> x(2);
    
    double minf = 0.0;

    x[0] = 0.25*M_PI;
    x[1] = 0.75*M_PI;
    
    opt.set_min_objective(myfunc, &parameters[0]);
    opt.optimize(x, minf);
   
    if (x[0] > M_PI)
        x[0] = 2.0*M_PI - x[0];
    if (x[1] > M_PI)
        x[1] = 2.0*M_PI - x[1];
    cout << "theta0 = ";
    cout << x[0]*180.0/M_PI << endl;
    cout << "theta1 = ";
    cout << x[1]*180.0/M_PI << endl;
    cout << "total energy = ";
    cout << minf << endl;
    */
    
    for(double field = 0.0; field<5.1;field+=0.1)
    {
        nlopt::opt opt(nlopt::LN_SBPLX,2);
        std::vector<double> ub(2);
        ub[0] = 2.0*M_PI;
        ub[1] = 2.0*M_PI;
        opt.set_upper_bounds(ub);
        std::vector<double> lb(2);
        lb[0] = 0.0;
        lb[0] = 0.0;
        opt.set_lower_bounds(lb);
        
        opt.set_ftol_rel(1.0e-10);
            
        std::vector<double> x(2);
        
        double minf = 0.0;
            
        x[0] = 0.75*M_PI;
        x[1] = 0.75*M_PI;
        
        double parameters[6] = {-2.5,-12.0,-12.0,0.0,-1.2,field};
        opt.set_min_objective(myfunc, &parameters[0]);
        opt.optimize(x, minf);
            
        cout << "H = " << field/2.0/5.7883818066e-2;
        cout << " Tesla";
        cout << ", Energy = " << minf;
        cout << " meV";
        cout << ", theta0 = ";
        if (x[0] > M_PI)
            x[0] = 2.0*M_PI - x[0];
        if (x[1] > M_PI)
            x[1] = 2.0*M_PI - x[1];
        cout << x[0]*180.0/M_PI;
        cout << " Degrees, theta1 = ";
        cout << x[1]*180.0/M_PI << " Degrees" << endl;
        
    }
 
    return 0;
}
