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
#include "Cell/Neighbors.h"
#include "Interactions/ExchangeInteraction.h"
#include "Interactions/AnisotropyInteraction.h"
#include "CalculateAngles.h"


using namespace std;
using namespace Eigen;

double SA = 2.1;
double SB = 0.7;

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
    } catch (...)
    {
        cout << "exception occurred" << endl;
    }
        
    return thetaphi[0];
}



double max_function(double sum, point x)
{
    if (x.frequency > 15.0 && x.frequency < 20.0)
        sum += x.intensity;
    return sum;
}

double min_function(double sum, point x)
{
    if (x.frequency < 15.0)
        sum += x.intensity;
    return sum;
}

double minimum(SpinWave &SW, double min_value, double max_value)
{
    double KX = -2.0;
    double KY = -2.0;
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
    
    double low_intensity = std::accumulate(points1.begin()+2,points1.end(),0.0,min_function);
    low_intensity += std::accumulate(points2.begin()+2,points2.end(),0.0,min_function);
    low_intensity += std::accumulate(points3.begin()+2,points3.end(),0.0,min_function);

    //cout << "low_intensity " << low_intensity << endl;
    
    double diff_sq = pow(low_intensity*1000.0,2);

    return diff_sq;
}

double maximum(SpinWave &SW, double min_value, double max_value)
{
    double KX = -2.0;
    double KY = -2.0;
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
    
    double high_intensity = std::accumulate(points1.begin()+2,points1.end(),0.0,max_function);
    high_intensity += std::accumulate(points2.begin()+2,points2.end(),0.0,max_function);
    high_intensity += std::accumulate(points3.begin()+2,points3.end(),0.0,max_function);
    
    //cout << "high_intensity " << high_intensity << endl;
    
    double diff_sq;
    
    diff_sq = 0.1/high_intensity;
    if (!isnormal(diff_sq))
    {
        diff_sq = 1000000000.0;
    }
    
    return diff_sq;
}

/*double max_function2(double sum, point x)
{
    if (x.frequency > 20.0 && x.frequency < 25.0)
        sum += x.intensity;
    return sum;
}

double min_function2(double sum, point x)
{
    if (x.frequency > 10.0 && x.frequency < 20.0)
        sum += x.intensity;
    return sum;
}

double minimum2(SpinWave &SW, double min_value, double max_value)
{
    double KX = -2.0;
    double KY = -2.0;
    double KZ = 1.0;
    
    SW.createMatrix(KX,KY,KZ);
    SW.Calc();
    vector<point> points1 = SW.getPoints();
    
    SW.createMatrix(KZ,KX,KY);
    SW.Calc();
    vector<point> points2 = SW.getPoints();
    
    SW.createMatrix(KY,KZ,KX);
    SW.Calc();
    vector<point> points3 = SW.getPoints();
    
    double low_intensity = std::accumulate(points1.begin(),points1.end(),0.0,min_function2);
    low_intensity += std::accumulate(points2.begin(),points2.end(),0.0,min_function2);
    low_intensity += std::accumulate(points3.begin(),points3.end(),0.0,min_function2);
    
    double diff_sq = low_intensity*1000.0;
    
    return diff_sq;
    
}

double maximum2(SpinWave &SW, double min_value, double max_value)
{
    double KX = -2.0;
    double KY = -2.0;
    double KZ = 1.0;
    
    SW.createMatrix(KX,KY,KZ);
    SW.Calc();
    vector<point> points1 = SW.getPoints();
    
    SW.createMatrix(KZ,KX,KY);
    SW.Calc();
    vector<point> points2 = SW.getPoints();
    
    SW.createMatrix(KY,KZ,KX);
    SW.Calc();
    vector<point> points3 = SW.getPoints();
    
    double high_intensity = std::accumulate(points1.begin(),points1.end(),0.0,max_function2);
    high_intensity += std::accumulate(points2.begin(),points2.end(),0.0,max_function2);
    high_intensity += std::accumulate(points3.begin(),points3.end(),0.0,max_function2);
    
    double diff_sq = 0.2/high_intensity;
    if (isnan(diff_sq))
    {
        diff_sq = 1000000000.0;
    }

    return diff_sq;
    
}
*/


double fitting(SpinWave &SW, double KX, double KY, double KZ, double freq, double error)
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
    
    
    //cout << "freq: " <<  freq << endl;
    /*for (auto it = points2.begin(); it != points2.end(); it++)
    {
        cout << it->frequency << " " << it->intensity << endl;
    }
    
    cout << endl;
     */
    int pos = 0;
    /*cout << points1[pos].frequency << " " << points1[pos].intensity << endl;
    cout << points1[1].frequency << " " << points1[1].intensity << endl;
    
    cout << points2[pos].frequency << " " << points2[pos].intensity << endl;
    cout << points2[1].frequency << " " << points2[1].intensity << endl;
    
    cout << points3[pos].frequency << " " << points3[pos].intensity << endl;
    cout << points3[1].frequency << " " << points3[1].intensity << endl;
    */
    double tmp0 = (points1[0].frequency+points2[0].frequency+points3[0].frequency)/3.0;
    double tmp1 = (points1[1].frequency+points2[1].frequency+points3[1].frequency)/3.0;
    
    if (pow(tmp0-freq,2)<pow(tmp1-freq,2))
        pos = 0;
    else
        pos = 1;

    double average = (points1[pos].frequency*points1[pos].intensity + points2[pos].frequency*points2[pos].intensity  + points3[pos].frequency*points3[pos].intensity);
    average = average/(points1[pos].intensity+points2[pos].intensity+points3[pos].intensity);
    //double average = points1[pos].frequency;
    //cout << average << endl;
    //cout << " " << endl;

    diff_sq += pow((average-freq)/error,2);
    
    return diff_sq;
}

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    double jbb = *(double *)my_func_data;
    //cout << "min_value = " << min_value << endl;
    /*cout << "x[i] = ";
    for (size_t i = 0; i<x.size();i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;
    */
    
    
    double theta = getCantingAngle(x[0],jbb,x[1]-jbb,x[2]);
    
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
    
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn0","V0",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn0","V1",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn0","V2",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn0","V3",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn1","V0",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn1","V1",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn1","V2",3.48,3.57));
    builder.Add_Interaction(new ExchangeInteraction("Jab",x[0],"Mn1","V3",3.48,3.57));
    
    builder.Add_Interaction(new ExchangeInteraction("Jbb",jbb,"V0","V1",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbb",jbb,"V2","V3",2.975,3.06));

    double jbbp = -1.0*jbb + x[1];
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",jbbp,"V0","V2",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",jbbp,"V0","V3",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",jbbp,"V1","V2",2.975,3.06));
    builder.Add_Interaction(new ExchangeInteraction("Jbbp",jbbp,"V1","V3",2.975,3.06));
    
    Vector3 direction(-1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[2],direction,"V0"));
    direction = Vector3(1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[2],direction,"V1"));
    direction = Vector3(1.0,1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dby",x[2],direction,"V2"));
    direction = Vector3(-1.0,-1.0,-1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Dbx",x[2],direction,"V3"));
    
    Vector3 zhat(0.0,0.0,1.0);
    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn0"));
    builder.Add_Interaction(new AnisotropyInteraction("Daz",x[3],zhat,"Mn1"));
    
    SpinWave SW = builder.Create_Element();
    
    double diff_sq = 0.0;
    
    /*diff_sq += fitting(SW,2.0,0.0,0.0898357,9.67548,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.300308,9.59135,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.490246,9.59135,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.703285,9.63341,1.0);
    diff_sq += fitting(SW,2.0,0.0,0.898357,9.54928,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.10370,8.96034,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.29877,7.86659,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.48871,6.35216,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.68121,4.33293,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.85062,2.52404,1.0);
    diff_sq += fitting(SW,2.0,0.0,1.94559,1.93510,0.5);
    diff_sq += fitting(SW,2.0,0.0,2.14836,2.52404,1.0);
    diff_sq += fitting(SW,2.0,0.0,2.31519,4.16466,1.0);
    diff_sq += fitting(SW,2.0,0.0,2.49743,6.05769,1.0);
    
    diff_sq += fitting(SW,4.0-1.09151,0.0,1.09151,9.16866,1.0);
    diff_sq += fitting(SW,4.0-1.29285,0.0,1.29285,7.78708,1.0);
    diff_sq += fitting(SW,4.0-1.49251,0.0,1.49251,6.57297,1.0);
    diff_sq += fitting(SW,4.0-1.65225,0.0,1.65225,5.35885,1.0);
    diff_sq += fitting(SW,4.0-1.75541,0.0,1.75541,4.31220,1.0);
    diff_sq += fitting(SW,4.0-1.84859,0.0,1.84859,3.01435,1.0);
    diff_sq += fitting(SW,4.0-1.94842,0.0,1.94842,1.92584,0.5);
    */
    
    /*diff_sq += fitting(SW,-2.0,-2.0,2.0,2.25,0.5);
    diff_sq += fitting(SW,-2.0,-2.0,1.75,3.69,0.25);
    diff_sq += fitting(SW,-2.0,-2.0,1.50,7.24,0.25);
    diff_sq += fitting(SW,-2.0,-2.0,1.25,9.19,0.25);
    diff_sq += fitting(SW,-2.0,-2.0,1.0,10.0,0.25);

    diff_sq += fitting(SW,-2.00,-2.00,0.0,2.15,0.5);
    diff_sq += fitting(SW,-2.25,-2.25,0.0,4.84,0.11);
    diff_sq += fitting(SW,-2.50,-2.50,0.0,7.48,0.15);
    diff_sq += fitting(SW,-2.75,-2.75,0.0,8.67,0.2);
    diff_sq += fitting(SW,-3.00,-3.00,0.0,10.0,0.25);
    
    diff_sq += fitting(SW,-1.0,-1.0,1.0,2.07,0.3);
    diff_sq += fitting(SW,-1.1,-1.1,0.9,3.77,0.18);
    diff_sq += fitting(SW,-1.2,-1.2,0.8,5.73,0.24);
    diff_sq += fitting(SW,-1.3,-1.3,0.7,7.79,0.17);
    diff_sq += fitting(SW,-1.4,-1.4,0.6,8.12,0.17);
    diff_sq += fitting(SW,-1.5,-1.5,0.5,8.53,0.3);
    */
    diff_sq += fitting(SW,-2.0,-2.0,2.0,2.25,1.0);
    diff_sq += fitting(SW,-2.0,-2.0,1.75,3.69,1.0);
    diff_sq += fitting(SW,-2.0,-2.0,1.50,7.24,1.0);
    diff_sq += fitting(SW,-2.0,-2.0,1.25,9.19,1.0);
    diff_sq += fitting(SW,-2.0,-2.0,1.0,10.0,1.0);
    
    diff_sq += fitting(SW,-2.00,-2.00,0.0,2.15,1.0);
    diff_sq += fitting(SW,-2.25,-2.25,0.0,4.84,1.0);
    diff_sq += fitting(SW,-2.50,-2.50,0.0,7.48,1.0);
    diff_sq += fitting(SW,-2.75,-2.75,0.0,8.67,1.0);
    diff_sq += fitting(SW,-3.00,-3.00,0.0,10.0,1.0);
    
    diff_sq += fitting(SW,-1.0,-1.0,1.0,2.07,1.0);
    diff_sq += fitting(SW,-1.1,-1.1,0.9,3.77,1.0);
    diff_sq += fitting(SW,-1.2,-1.2,0.8,5.73,1.0);
    diff_sq += fitting(SW,-1.3,-1.3,0.7,7.79,1.0);
    diff_sq += fitting(SW,-1.4,-1.4,0.6,8.12,1.0);
    diff_sq += fitting(SW,-1.5,-1.5,0.5,8.53,1.0);
    
    
    /*diff_sq += fitting(SW,-2.32362,-2.32362,0.0,6.11621,1.0);
    diff_sq += fitting(SW,-2.26214,-2.26214,0.0,5.12232,1.0);
    diff_sq += fitting(SW,-2.18770,-2.18770,0.0,4.05199,1.0);
    diff_sq += fitting(SW,-2.10032,-2.10032,0.0,2.79052,1.0);
    diff_sq += fitting(SW,-1.97735,-1.97735,0.0,1.83486,1.0);
    diff_sq += fitting(SW,-1.88997,-1.88997,0.0,2.98165,1.0);
    diff_sq += fitting(SW,-1.80259,-1.80259,0.0,4.09021,1.0);
    diff_sq += fitting(SW,-1.74757,-1.74757,0.0,5.12232,1.0);
    diff_sq += fitting(SW,-1.64078,-1.64078,0.0,6.23089,1.0);
    diff_sq += fitting(SW,-1.57282,-1.57282,0.0,7.30122,1.0);
    diff_sq += fitting(SW,-1.46926,-1.46926,0.0,8.33333,1.0);
    diff_sq += fitting(SW,-1.28803,-1.28803,0.0,8.67737,1.0);
    diff_sq += fitting(SW,-1.20065,-1.20065,0.0,9.51835,1.0);
    diff_sq += fitting(SW,-1.07767,-1.07767,0.0,9.63303,1.0);
    diff_sq += fitting(SW,-0.996764,-0.996764,0.0,9.86238,1.0);
    diff_sq += fitting(SW,-0.877023,-0.877023,0.0,9.78593,1.0);
    
    diff_sq += fitting(SW,0.0,0.0,-0.795729,9.62897,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.699164,9.05477,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.595172,8.65724,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.502321,7.64134,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.424327,6.80212,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.379759,6.31625,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.309192,5.38869,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.257196,4.41696,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.190344,3.84276,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.0900650,3.09187,1.0);
    diff_sq += fitting(SW,0.0,0.0,-0.0157845,2.56184,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.0547817,2.34099,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.129062,3.18021,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.173630,3.71025,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.207057,4.10777,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.292479,5.03534,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.351903,5.69788,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.429898,6.53710,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.507892,7.15548,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.604457,8.03887,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.671309,8.83392,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.771588,9.45230,1.0);
    diff_sq += fitting(SW,0.0,0.0,0.875580,9.49647,1.0);
    diff_sq += fitting(SW,0.0,0.0,1.00928,9.49647,1.0);
*/
    
    
    //cout << "diff_sq before = " << diff_sq << endl;
    diff_sq += minimum(SW,15.0,20.0);
    diff_sq += maximum(SW,15.0,20.0);
    //diff_sq += minimum2(SW,20.0,25.0);
    //diff_sq += maximum2(SW,20.0,25.0);


    double CantingAngle = 180.0/M_PI*(M_PI - theta);
    //cout << "Canting Angle= " << CantingAngle << endl;
    diff_sq += 1.0*pow(CantingAngle-36.0,2);
    //cout << "diff_sq after = " << diff_sq << endl;
    
    //cout << "diff_sq= " << diff_sq << endl;
    return diff_sq;
}

int main()
{
    cout << "#min_value JAB JBB JBBP DB111 SA/SB DAZ error" << endl;

    std::vector<double> ub(4);
    ub[0] = -1.0;
    ub[1] = -4.0;
    ub[2] = -7.0;
    ub[3] =  1.0;
    
    std::vector<double> lb(4);
    lb[0] = -3.0;
    lb[1] = -10.0;
    lb[2] = -12.0;
    lb[3] =  -1.0;

    std::vector<double> x(4);
    
    double minf = 0.0;
    
    //vector<double> gradient;

    for(double min_value = -30.0; min_value < -9.1 ;min_value += 1.0)
    {

        x[0] = -1.8;
        x[1] = -6.9;
        x[2] = -9.1;
        x[3] =  0.4;

        nlopt::opt opt(nlopt::LN_SBPLX,4);
        opt.set_min_objective(myfunc, &min_value);
        opt.set_upper_bounds(ub);
        opt.set_lower_bounds(lb);
        opt.set_ftol_abs(5.0e-6);
        opt.set_maxeval(100000);
        opt.optimize(x, minf);
    
        //cout << result << endl;
    
        //cout << minf << endl;
        
        cout << min_value << " ";
        for (vector<double>::iterator it = x.begin(); it!=x.end();it++)
            cout << (*it) << " ";
        
        double theta = getCantingAngle(x[0],min_value,x[1]-min_value,x[2]);
        cout << 180.0/M_PI*(M_PI-theta) << " ";
        cout << minf << endl;

    }

    return 0;
}
