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

int main(int argc, const char * argv[])
{
    vector<double> X(4,0.0);
    vector<SW_sublattice> SL(2);
    double KX;
    SW_Matrix test;
    
    Vector3d pos;
    vector<Vector3d> positions;
    vector<int> types;

    types.push_back(1);
    pos << 0.0,0.0,0.5;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.5,-0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
    pos << 0.0,0.0,-0.5;
    positions.push_back(pos);
    types.push_back(1);
    pos << -0.5,0.5,0.0;
    positions.push_back(pos);
    types.push_back(1);
    pos << -0.5,-0.5,0.0;
    positions.push_back(pos);
    
    SL[0].set_sublattice(5.0/2.0,M_PI/2.0 - 0.01098,0.0);
    SL[0].add_neighbors(types,positions);
    
    vector<double> interaction;
    interaction.push_back(0.0);
    interaction.push_back(1.0/sqrt(2.0));
    SL[0].set_interactions(interaction);
    
    SL[1].set_sublattice(5.0/2.0,M_PI/2.0 - 0.01098,M_PI);
    for (int i=0;i<6;i++)
    {
        types[i] = 0;
    }
    
    SL[1].add_neighbors(types,positions);
    
    interaction.erase(interaction.begin(),interaction.end());
    
    interaction.push_back(1.0/sqrt(2.0));
    interaction.push_back(0.0);
    
    SL[1].set_interactions(interaction);
    
    //cout << SL[0].get_rot_matrix() << endl;
    //cout << SL[1].get_rot_matrix() << endl;
    
    X[0] = 2.112; //2.48;
    X[1] = 0.109;
    X[2] = 4.61e-3;
    X[3] = 1.135e-3;
        
    test.set_parr(SL,X);
    for (int i=0;i<301;i++)
    {
        KX = (double)i/100.0 - 1.5;
        //cout << KX << endl;
        test.CreateMatrix_YFeO3(0.0,KX,-3.0);
        test.Calc_Eigenvalues();
        test.Calc_Weights();
        test.Rotation_Matrix();
        test.Unique_Solutions();
        test.Signif_Solutions(0.0,KX,-3.0);
        //cout << endl;
    }
    return 0;
}
