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
#include "Interactions/Exch_Interaction.h"
#include "Interactions/AnisotropyInteraction.h"


using namespace std;
using namespace Eigen;

double minimum(SpinWave &SW, double min_value)
{
    double KX = 2.0;
    double KY = 0.0;
    double KZ = 2.0;
    
    SW.createMatrix(KX,KY,KZ);
    SW.Calc();
    vector<point> points1 = SW.getPoints();
    
    SW.createMatrix(KZ,KX,KY);
    SW.Calc();
    vector<point> points2 = SW.getPoints();
    
    SW.createMatrix(KY,KZ,KX);
    SW.Calc();
    vector<point> points3 = SW.getPoints();
    
    /*cout << points1[2].frequency << "  " << points2[2].frequency << "  " << points3[2].frequency << endl;
    cout << points1[2].intensity << "  " << points2[2].intensity << "  " << points3[2].intensity << endl;
    cout << points1[3].frequency << "  " << points2[3].frequency << "  " << points3[3].frequency << endl;
    cout << points1[3].intensity << "  " << points2[3].intensity << "  " << points3[3].intensity << endl;
    cout << points1[4].frequency << "  " << points2[4].frequency << "  " << points3[4].frequency << endl;
    cout << points1[4].intensity << "  " << points2[4].intensity << "  " << points3[4].intensity << endl;
    cout << points1[5].frequency << "  " << points2[5].frequency << "  " << points3[5].frequency << endl;
    cout << points1[5].intensity << "  " << points2[5].intensity << "  " << points3[5].intensity << endl;
    cout << endl;
    */
    double lower_limit = (points1[2].frequency*points1[2].intensity + points2[2].frequency*points2[2].intensity  + points3[2].frequency*points3[2].intensity)/(points1[2].intensity+points2[2].intensity+points3[2].intensity);

    double diff_sq = pow(lower_limit-min_value,2);
    return diff_sq;

}


double fitting(SpinWave &SW, double KX, double KY, double KZ, double freq, long pos, double error)
{
    
    double diff_sq = 0.0;
    SW.createMatrix(KX,KY,KZ);
    SW.Calc();
    vector<point> points1 = SW.getPoints();
    
    SW.createMatrix(KZ,KX,KY);
    SW.Calc();
    vector<point> points2 = SW.getPoints();
    
    SW.createMatrix(KY,KZ,KX);
    SW.Calc();
    vector<point> points3 = SW.getPoints();
    
    /*cout << "freq: " <<  freq << endl;
    
    cout << points1[0].frequency << endl;
    cout << points1[1].frequency << endl;
    
    cout << points2[0].frequency << endl;
    cout << points2[1].frequency << endl;
    
    cout << points3[0].frequency << endl;
    cout << points3[1].frequency << endl;

    cout << " " << endl;
    */
    double average = (points1[pos].frequency*points1[pos].intensity + points2[pos].frequency*points2[pos].intensity  + points3[pos].frequency*points3[pos].intensity);
    average = average/(points1[pos].intensity+points2[pos].intensity+points3[pos].intensity);
    diff_sq += pow((average-freq)/error,2);
    
    return diff_sq;
}

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    double min_value = *(double *)my_func_data;

    /*cout << "x[i] = ";
    for (size_t i = 0; i<6;i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    */
    
    double SA = 1.54;
    double SB = 0.9;
    double theta = M_PI - 35.0*M_PI/180.0;
    
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

    builder.Add_Interaction(new Exch_Interaction("Jaa",x[0],"Mn0","Mn1",3.0,5.0));
    
    builder.Add_Interaction(new Exch_Interaction("Jbb",x[1],"V0","V1",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbb",x[1],"V2","V3",2.975,3.06));

    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V0","V2",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V0","V3",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V1","V2",2.975,3.06));
    builder.Add_Interaction(new Exch_Interaction("Jbbp",x[2],"V1","V3",2.975,3.06));
    
    Vector3 zhat(0.0,0.0,1.0);

    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn0"));
    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn1"));
    
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V0"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V1"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V2"));
    builder.Add_Interaction(new AnisotropyInteraction("Dbz",x[4],zhat,"V3"));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[5],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[5],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[5],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[5],direction,"V3"));

    
    double JBB = x[1];
    double JBBP = x[2];
    double DB = x[5];
    double DBz = x[4];
    double JAB = SB*((6.0*JBB+6.0*JBBP+DB-3.0*DBz)*cos(theta)*sin(theta) -sqrt(2.0)*DB*(2.0*pow(cos(theta),2)-1.0))/(-9.0*SA*sin(theta));
    
    //cout << "JAB= " << JAB << endl;
    
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V0",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V1",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V2",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn0","V3",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V0",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V1",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V2",3.48,3.57));
    builder.Add_Interaction(new Exch_Interaction("Jab",JAB,"Mn1","V3",3.48,3.57));
    
    SpinWave SW = builder.Create_Element();
    
    double diff_sq = 0.0;
    
    /*//MnV2O4
    diff_sq += fitting(SW,2.0,0.0,0.0,9.67548,1,1.0);
    diff_sq += 2.0*fitting(SW,2.0,0.0,2.0,2.0,0,1.0);
    diff_sq += fitting(SW,4.0-1.0,0.0,1.0,9.16866,0,1.0);
    diff_sq += minimum(SW, min_value);
    */
    //Mn0.6Co0.4V2O4
    ///diff_sq += fitting(SW,1.0,1.0,0.0,10.0,0,1.0);
    ///diff_sq += fitting(SW,2.0,2.0,0.0,2.0,0,0.1);
    //diff_sq += minimum(SW, min_value);
 
    /*//Mn0.4=6Co0.4V2O4
    diff_sq += fitting(SW,0.0,0.0,4.00,2.0,0,1.0);
    //diff_sq += fitting(SW,0.0,0.0,3.75,5.7,0,1.0,min_value);
    //diff_sq += fitting(SW,0.0,0.0,3.50,9.23,0,1.0,min_value);
    //diff_sq += fitting(SW,0.0,0.0,3.25,9.92,0,1.0,min_value);
    diff_sq += fitting(SW,0.0,0.0,3.00,10.66,0,1.0);

    diff_sq += fitting(SW,2.00,2.00,0.0,2.0,0,0.5);
    //diff_sq += fitting(SW,1.75,1.75,0.0,6.96,0,1.0,min_value);
    //diff_sq += fitting(SW,1.50,1.50,0.0,8.12,0,1.0,min_value);
    //diff_sq += fitting(SW,1.25,1.25,0.0,9.73,0,1.0,min_value);
    diff_sq += fitting(SW,1.00,1.00,0.0,10.49,0,1.0);
    
    //diff_sq += fitting(SW,1.85,1.85,-1.85,4.87,0,1.0);
    //diff_sq += fitting(SW,1.75,1.75,-1.75,7.54,0,1.0);
    //diff_sq += fitting(SW,1.50,1.50,-1.50,8.3,0,1.0);
    //diff_sq += fitting(SW,1.00,1.00,-1.00,8.33,0,1.0);
*/
    
    diff_sq += fitting(SW,2.0,2.0,0.0,2.0,0,1.0);
    //cout << diff_sq << endl;
    diff_sq += fitting(SW,3.0,3.0,0.0,11.5,1,1.0);
    //cout << diff_sq << endl;
    diff_sq += minimum(SW, min_value);
    //cout << diff_sq << endl;

    
    //cout << "diff_sq= " << diff_sq << endl;
    return diff_sq;
}

int main()
{
    nlopt::opt opt(nlopt::LN_COBYLA,6);
    std::vector<double> ub(6);
    ub[0] =  0.0;
    ub[1] =-10.0;
    ub[2] = 10.0;
    ub[3] =  2.0;
    ub[4] =  0.0;
    ub[5] = -2.0;
    opt.set_upper_bounds(ub);
    std::vector<double> lb(6);
    lb[0] =  0.0;
    lb[1] = -10.0;
    lb[2] =  5.0;
    lb[3] =  0.3;
    lb[4] =  0.0;
    lb[5] = -12.0;
    opt.set_lower_bounds(lb);

    opt.set_stopval(1.0e-2);
    opt.set_maxeval(500);

    std::vector<double> x(6);
    
    double minf = 0.0;
    
    cout << "#MinValue JAB JBB JBBP DAZ DBZ DB111 " << endl;
    
    vector<double> gradient;
    /*
    x[0] = 0.0;
    x[1] =-12.0;
    x[2] = 5.0;
    x[3] = 0.2;
    x[4] = 0.0;
    x[5] = -5.0;
    */
    
    x[0] = 0.0;
    x[1] =-10.0;
    //x[2] =  5.0;
    x[2] = 5.0;
    //x[3] = 1.0;
    x[3] = 0.9;
    x[4] = 0.0;
    x[5] = -4.5;
    
    for(double min_value = 14.0; min_value < 21.1 ;min_value += 0.1)
    {


        
        opt.set_min_objective(myfunc, &min_value);
        opt.optimize(x, minf);
    
        //cout << result << endl;
    
        //cout << minf << endl;
        
        double SA = 1.54;
        double SB = 0.9;
        double theta = M_PI - 35.0*M_PI/180.0;
        double JBB = x[1];
        double JBBP = x[2];
        double DB = x[5];
        double DBz = x[4];
        double JAB = SB*((6.0*JBB+6.0*JBBP+DB-3.0*DBz)*cos(theta)*sin(theta) -sqrt(2.0)*DB*(2.0*pow(cos(theta),2)-1.0))/(-9.0*SA*sin(theta));
        
        cout << min_value << " " << JAB << " ";
        for (vector<double>::iterator it = x.begin()+1; it!=x.end();it++)
            cout << (*it) << " ";
        
        cout << myfunc(x,gradient,&min_value) << endl;
        //cout << endl;
        //cout << "chi_squared = " << minf << endl;
        //cout << "reduced chi_squared = " << minf/21.0 << endl;
    }
    /*vector<double> gradient;
    

    //cout << endl;
    long i = 4;
    vector<double> xtemp = x;
    for (int j=-20;j!=21;j++)
    {
        xtemp[i] = x[i] + x[i]*(double)j/0.2;
        double diffsq = myfunc(xtemp,gradient,NULL);
        cout << "x[" << i << "] = " << xtemp[i];
        cout << ", reduced chi^2= " << (diffsq - minf)/21.0 -0.0 << endl;
    }*/

    return 0;
}
