//
//  main.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "SW_Matrix.h"
#include "SW_sublattice.h"

using namespace Eigen;
using namespace std;


struct Results
{
    vector<double> kpoint;
    vector<double> frequencies;
    vector<double> intensities;
};

int main(int argc, const char * argv[])
{
    vector<double> X(6,0.0);
    vector<SW_sublattice> SL(4);
    double KX,KY,KZ;
    SW_Matrix test;
    vector<Results> SW_results;
    
    Vector3d pos;
    vector<Vector3d> positions;
    vector<int> types;
    
    
    X[0] = -4.984; //-2.48*2.0;
    X[1] = 0.0; //0.109*2.0;
    X[2] = 1.135e-3;
    X[3] = 4.61e-3;
    X[4] = 0.109;
    X[5] = -0.109;
    
    double delta = 2.0*X[5]/(6*X[0]+X[3]-X[2]);
    
    //double delta = M_PI/12.0;
    double phi = 0.0;
    cout << "delta= " << delta << endl;
    
    //nn
    types.push_back(1);
    pos << 0.0,0.0,0.5;
    positions.push_back(pos);
    types.push_back(3);
    pos << 0.5,-0.5,0.0;
    positions.push_back(pos);
    types.push_back(3);
    pos << 0.5,0.5,0.0;
    positions.push_back(pos);
    //nnn
    types.push_back(0);
    pos << 1.0,0.0,0.0;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.0,1.0,0.0;
    positions.push_back(pos);
    types.push_back(2);
    pos << 0.5,0.5,0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << 0.5,-0.5,0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << -0.5,0.5,0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << -0.5,-0.5,0.5;
    positions.push_back(pos);

    SL[0].set_sublattice(5.0/2.0,M_PI/2.0 - delta,M_PI);
    SL[0].add_neighbors(types,positions);
    
    SL[1].set_sublattice(5.0/2.0,M_PI/2.0 - delta,0.0);
    for (int i=0;i<1;i++)
    {
        types[i] = 0;
    }
    for (int i=1;i<3;i++)
    {
        types[i] = 2;
    }
    for (int i=3;i<5;i++)
    {
        types[i] = 1;
    }
    for (int i=5;i<9;i++)
    {
        types[i] = 3;
    }
    
    //for (int i=0;i<18;i++)
    //    cout << types[i] << endl;
    
    SL[1].add_neighbors(types,positions);
    
    SL[2].set_sublattice(5.0/2.0,M_PI/2.0 - delta,M_PI);
    for (int i=0;i<1;i++)
    {
        types[i] = 3;
    }
    for (int i=1;i<3;i++)
    {
        types[i] = 1;
    }
    for (int i=3;i<5;i++)
    {
        types[i] = 2;
    }
    for (int i=5;i<9;i++)
    {
        types[i] = 0;
    }
    
    SL[2].add_neighbors(types,positions);
    
    SL[3].set_sublattice(5.0/2.0,M_PI/2.0 - delta,0.0);
    for (int i=0;i<1;i++)
    {
        types[i] = 2;
    }
    for (int i=1;i<3;i++)
    {
        types[i] = 0;
    }
    for (int i=3;i<5;i++)
    {
        types[i] = 3;
    }
    for (int i=5;i<9;i++)
    {
        types[i] = 1;
    }
    
    SL[3].add_neighbors(types,positions);
    
    cout << SL[0].get_rot_matrix() << endl;
    cout << SL[1].get_rot_matrix() << endl;
    cout << SL[2].get_rot_matrix() << endl;
    cout << SL[3].get_rot_matrix() << endl;

    KX = 0.0;
    KZ = 0.0;
    //KX = 1.0;
    //KZ = -3.0;
    vector<double> kpoint (3);
    Results tmp_result;
    for (int i=0;i<301;i++)
    {
        KY = (double)i/100.0 - 1.5;
        cout << KY << endl;
        test.set_parr(SL,X);
	    test.CreateMatrix_exchange(KX,KY,KZ);
        test.CreateMatrix_anis_z();
        test.CreateMatrix_anis_x();
        test.CreateMatrix_DMy(KX,KY,KZ);
        test.Calc_Eigenvalues();
        test.Calc_Weights();
        test.Rotation_Matrix();
        test.Unique_Solutions();
        test.Signif_Solutions(KX,KY,KZ);
        kpoint[0] = KX;
        kpoint[1] = KY;
        kpoint[2] = KZ;
        tmp_result.kpoint = kpoint;
        tmp_result.frequencies = test.Get_Frequencies();
        tmp_result.intensities = test.Get_Intensities();
        SW_results.push_back(tmp_result);
    }
    
    for (int i=0;i<301;i++)
    {
        cout << SW_results[i].kpoint[0] << '\t' << SW_results[i].kpoint[1] << '\t' << SW_results[i].kpoint[2] << '\t';
        for (int k=0; k!=SW_results[i].frequencies.size();k++)
        {
            cout << SW_results[i].frequencies[k] << '\t';
            
        }
        for (int k=0; k!=SW_results[i].intensities.size();k++)
        {
            cout << SW_results[i].intensities[k] << '\t';
            
        }
        
        cout << sqrt(12.0*abs(X[0])*6.25*(2.0*(X[3]-X[2]))) << '\t';
        cout << sqrt(12.0*abs(X[0])*6.25*(6.0*abs(X[5])*tan(delta) + 2.0*X[3])) << '\t';
        cout << sqrt(12.0*abs(X[0])*6.25*(-6.0/3.0*abs(X[5])*tan(delta) + 2.0*X[3])) << '\t';
        cout << endl;
        
        
        
    }
    return 0;
}
