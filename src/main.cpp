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
    double KX,KY,KZ;
    SW_Matrix test;
    
    X[0] = 1.0; //J
    X[1] = -2.0; //g*J
    X[2] = 2.0; //-n*J
    X[3] = -8.0; //B'
    
    KX = 0.0;
    KZ = 0.0;
    
    Vector3d pos;
    vector<Vector3d> positions;
    vector<int> types;

    double alpha = acos(X[3]/(4.0*X[1]/X[0]));
    cout << "alpha= " << alpha <<endl;
    
    types.push_back(1);
    pos << 0.0,0.5,0.0;
    positions.push_back(pos);
    //types.push_back(1);
    //pos << 0.0,-0.5,0.0;
    //positions.push_back(pos);
    types.push_back(0);
    pos << 1.0,0.0,0.0;
    positions.push_back(pos);
    //types.push_back(0);
    //pos << -1.0,0.0,0.0;
    //positions.push_back(pos);

    SL[0].set_sublattice(1.0,alpha,0.0);
    SL[0].add_neighbors(types,positions);
    
    SL[1].set_sublattice(1.0,alpha,M_PI);
    for (int i=0;i<1;i++)
    {
        types[i] = 0;
    }
    for (int i=1;i<2;i++)
    {
    	types[i] = 1;
    }
    SL[1].add_neighbors(types,positions);

    for (int i=0;i<11;i++)
    {
        KY = (double)i/10.0;
        //cout << KX << endl;
        test.set_parr(SL,X);
        test.CreateMatrix_exchange(KX,KY,KZ);
        test.CreateMatrix_bfield();
        //test.CreateMatrix_anis_z();
        //test.CreateMatrix_anis_x();
        //test.CreateMatrix_AFM(0.0,0.0,0.0);
        test.Calc_Eigenvalues();
        test.Calc_Weights();
        test.Rotation_Matrix();
        test.Unique_Solutions();
        test.Signif_Solutions(KX,KY,KZ);
        
        KY=KY/2.0;
        double S=SL[0].get_sublattice()[0];
        double J = X[0];
        double gamma = X[1]/X[0];
        double eta = -1.0*X[2]/X[0];
        double R2p = (eta-1.0)*(cos(2.0*M_PI*KX) -1.0) + 2.0*gamma*(-1.0+cos(2.0*M_PI*KY));
        double R2m = (eta-1.0)*(cos(2.0*M_PI*KX) -1.0) + 2.0*gamma*(-1.0-cos(2.0*M_PI*KY));
        
        cout << "analytical results" << endl;
        cout << KY << '\t' << sqrt(R2p*R2m)+(eta+1)*(1.0-cos(2.0*M_PI*KX));
        cout << '\t' << sqrt(R2p/R2m) << endl;
        cout << KY << '\t' << sqrt(R2p*R2m)-(eta+1)*(1.0-cos(2.0*M_PI*KX));
        cout << '\t' << sqrt(R2p/R2m) << endl;                     
        //cout << endl;
    }
    return 0;
}
