//
//  main.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
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
    for (int a=0;a<1;a++)
    {
        
    vector<double> X(5,0.0);
    vector<SW_sublattice> SL(2);
    double KX,KY,KZ;
    SW_Matrix test;
    vector<Results> SW_results;
        
    Vector3d pos;
    vector<Vector3d> positions;
    vector<int> types;
    
    X[0] = -2.505*2.0;
    X[1] = -0.162*2.0;
    X[2] = 3.47e-3;
    X[3] = 5.97e-3;
    X[4] = 0.0;
        
    //X[0] = -5.0;
    //X[1] = 0.0;
    //X[2] = 1.0e-3;
    //X[3] = 10.0e-3;

        
        
    
    //double delta= -2.0*X[4]/(6.0*X[0]+X[2]-X[3]);
    double delta = 0.01;
    
    X[4] = -0.5*delta*(6.0*X[0]+X[2]-X[3]);

    cout << "delta= " << delta << endl;
    cout << "DMy= " << X[4] << endl;

    //nn
    types.push_back(1);
    pos << 0.0,0.0,0.5;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.0,0.0,-0.5;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.5,-0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
    pos << -0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
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
    types.push_back(0);
    pos << 0.5,0.5,0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.5,-0.5,0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << -0.5,0.5,0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << -0.5,-0.5,0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << -0.5,-0.5,-0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << -0.5,0.5,-0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.5,-0.5,-0.5;
    positions.push_back(pos);
    types.push_back(0);
    pos << 0.5,0.5,-0.5;
    positions.push_back(pos);

    SL[0].set_sublattice(5.0/2.0,M_PI/2.0 - delta,M_PI);
    SL[0].add_neighbors(types,positions);
    
        SL[1].set_sublattice(5.0/2.0,M_PI/2.0 - delta,0.0);
    for (int i=0;i<6;i++)
    {
        types[i] = 0;
    }
    for (int i=6;i<18;i++)
    {
        types[i] = 1;
    }
    
    SL[1].add_neighbors(types,positions);
    
    MatrixXi exch_interaction, DMy_interaction, DMz_interaction;
    
    test.set_parr(SL,X);    
    exch_interaction.resize(2,2);
    exch_interaction << 1,0,
    0,1;
        
    test.set_exch_int(exch_interaction);
    
    DMy_interaction.resize(2,2);
    DMy_interaction << -1,4,
    -1,-1;
    
    test.set_DMy_int(DMy_interaction);
        
    //cout << SL[0].get_rot_matrix() << endl;
    //cout << SL[1].get_rot_matrix() << endl;
    //cout << SL[2].get_rot_matrix() << endl;
    //cout << SL[3].get_rot_matrix() << endl;

    KX = 1.0;
    KY = 0.0;
    KZ = -3.0;//(double)a;
        
    //for (int i=0;i<101;i++)
    //    cout << (double)rand()/RAND_MAX << endl;
        
    
    vector<double> kpoint (3);
    Results tmp_result;
            test.set_parr(SL,X);
    for (int i=0;i<151;i++)
    {
        KY = (double)i/50.0 + -1.5;
        cout << KX << '\t' << KY << '\t' << KZ << endl;
        test.Clear_LN();
        test.Calc_SW(KX,KY,KZ);
        kpoint[0] = KX;
        kpoint[1] = KY;
        kpoint[2] = KZ;
        tmp_result.kpoint = kpoint;
        tmp_result.frequencies = test.Get_Frequencies();
        tmp_result.intensities = test.Get_Intensities();
        SW_results.push_back(tmp_result);
    }
    
    ostringstream strs;
    strs << a;
    string output = "DMy_" + strs.str() + ".txt";
    ofstream out_file ( output.c_str() );
    
        for (int i=0;i<151;i++)
    {
        out_file << SW_results[i].kpoint[0] << '\t' << SW_results[i].kpoint[1] << '\t' << SW_results[i].kpoint[2] << '\t';
        for (int k=0; k!=SW_results[i].frequencies.size();k++)
        {
            out_file << SW_results[i].frequencies[k] << '\t';
            
        }
        for (int k=0; k!=SW_results[i].intensities.size();k++)
        {
            out_file << SW_results[i].intensities[k] << '\t';
            
        }
        
        out_file << endl;
    }
    out_file.close();
    }
    return 0;
}
