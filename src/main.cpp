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
    
    
    X[0] = -2.505*2.0;
    X[1] = -0.162*2.0;
    X[2] = 1.135e-3;
    X[3] = 4.61e-3;
    X[4] = 0.109;
    X[5] = 0.109;
    
    double delta= -2.0*X[4]/(6.0*X[0]+X[2]-X[3]);
    double phi = -2.0*X[5]/(4.0*X[0] - 8.0*X[1] - X[3]);
    //delta = 0.0;
    //phi = 0.0;
    cout << "delta= " << delta << endl;
    cout << "phi= " << phi << endl;
    
    
    //nn
    types.push_back(1);
    pos << 0.0,0.0,0.5;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.0,0.0,-0.5;
    positions.push_back(pos);
    types.push_back(3);
    pos << 0.5,-0.5,0.0;
    positions.push_back(pos);
    types.push_back(3);
    pos << 0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(3);
    pos << -0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(3);
    pos << -0.5,-0.5,0.0;
    positions.push_back(pos);
    //nnn
    types.push_back(0);
    pos << 1.0,0.0,0.0;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.0,1.0,0.0;
    positions.push_back(pos);
    types.push_back(0);
    pos << -1.0,0.0,0.0;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.0,-1.0,0.0;
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
    types.push_back(2);
    pos << -0.5,-0.5,-0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << -0.5,0.5,-0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << 0.5,-0.5,-0.5;
    positions.push_back(pos);
    types.push_back(2);
    pos << 0.5,0.5,-0.5;
    positions.push_back(pos);

    SL[0].set_sublattice(5.0/2.0,M_PI/2.0 - delta,M_PI+phi);
    SL[0].add_neighbors(types,positions);
    
    SL[1].set_sublattice(5.0/2.0,M_PI/2.0 - delta,0.0+phi);
    for (int i=0;i<2;i++)
    {
        types[i] = 0;
    }
    for (int i=2;i<6;i++)
    {
        types[i] = 2;
    }
    for (int i=6;i<10;i++)
    {
        types[i] = 1;
    }
    for (int i=10;i<18;i++)
    {
        types[i] = 3;
    }
    
    SL[1].add_neighbors(types,positions);
    
    SL[2].set_sublattice(5.0/2.0,M_PI/2.0 - delta,M_PI-phi);
    for (int i=0;i<2;i++)
    {
        types[i] = 3;
    }
    for (int i=2;i<6;i++)
    {
        types[i] = 1;
    }
    for (int i=6;i<10;i++)
    {
        types[i] = 2;
    }
    for (int i=10;i<18;i++)
    {
        types[i] = 0;
    }
    
    SL[2].add_neighbors(types,positions);
    
    SL[3].set_sublattice(5.0/2.0,M_PI/2.0 - delta,2.0*M_PI-phi);
    for (int i=0;i<2;i++)
    {
        types[i] = 2;
    }
    for (int i=2;i<6;i++)
    {
        types[i] = 0;
    }
    for (int i=6;i<10;i++)
    {
        types[i] = 3;
    }
    for (int i=10;i<18;i++)
    {
        types[i] = 1;
    }
    
    SL[3].add_neighbors(types,positions);
    
    //cout << SL[0].get_rot_matrix() << endl;
    //cout << SL[1].get_rot_matrix() << endl;
    //cout << SL[2].get_rot_matrix() << endl;
    //cout << SL[3].get_rot_matrix() << endl;

    KX = 1.0;
    KZ = -5.0;
    
    vector<double> kpoint (3);
    Results tmp_result;
    for (int i=0;i<151;i++)
    {
        //KX = 0.5 - (double)i/200.0 ;
        KY = -1.5 + (double)i/50.0;
        cout << KY << endl;
        test.set_parr(SL,X);
	    test.CreateMatrix_exchange_sumoverhalf(KX,KY,KZ);
        test.CreateMatrix_anis_z();
        test.CreateMatrix_anis_x();
        test.CreateMatrix_DMy(KX,KY,KZ);
        test.CreateMatrix_DMz(KX,KY,KZ);
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
    
    for (int i=0;i<151;i++)
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
        
        cout << endl;
    }
    return 0;
}
