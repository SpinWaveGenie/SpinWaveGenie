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
#include <iomanip>
#include <vector>
#include "SW_Matrix.h"
#include "SW_sublattice.h"

using namespace Eigen;
using namespace std;

int main(int argc, const char * argv[])
{
    vector<double> X(2,0.0);
    vector<SW_sublattice> SL(2);
    double KX,KY,KZ;
    SW_Matrix test;
    
    Vector3d pos;
    vector<Vector3d> positions;
    vector<int> types;

    types.push_back(1);
    pos << 1.0,0.0,0.0;
    positions.push_back(pos);
    
    SL[0].set_sublattice(1.0,0.0,0.0);
    SL[0].add_neighbors(types,positions);
    
    SL[1].set_sublattice(1.0,M_PI,0.0);
    for (int i=0;i<1;i++)
    {
        types[i] = 0;
    }
    SL[1].add_neighbors(types,positions);
    
    X[0] = -1.0;
    X[1] = -2.0;
        
    test.set_parr(SL,X);
    
    KY = 0.0;
    KZ = 0.0;
    
    for (int i=0;i<11;i++)
    {
        KX = (double)i/100.0;
        //cout << KX << endl;
        test.set_parr(SL,X);
	    test.CreateMatrix_exchange(KX,KY,KZ);
        //test.CreateMatrix_anis_z();
        test.CreateMatrix_anis_x();
        //test.CreateMatrix_AFM(0.0,0.0,0.0);
        test.Calc_Eigenvalues();
        test.Calc_Weights();
        test.Rotation_Matrix();
        test.Unique_Solutions();
        test.Signif_Solutions(KX,KY,KZ);
        
        double z=1.0;
        double S=SL[0].get_sublattice()[0];
        double J = abs(X[0]);
        double D = abs(X[1]);

        cout << "Analytical Solution" << endl;
        cout << scientific;
        cout << KX << '\t';
        cout << setprecision (9) << sqrt(4.0*pow(J*S,2)*(1.0-pow(cos(2.0*M_PI*KX),2))+4*J*D*pow(S,2)*(1+cos(2.0*M_PI*KX))) << '\t';
        cout << setprecision (9) << S/(4.0*sqrt(2.0*J*S))/sqrt(1.0+cos(2.0*M_PI*KX))*sqrt(2.0*J*S*(1.0-cos(2.0*M_PI*KX))+2.0*D*S) << endl;
        cout << 0.25*sqrt(D/2.0/J)*S << endl;
        cout << KX << '\t';
        cout << setprecision (9) << sqrt(4.0*pow(J*S,2)*(1.0-pow(cos(2.0*M_PI*KX),2))+4*J*D*pow(S,2)*(1-cos(2.0*M_PI*KX)))<< '\t';
        cout << setprecision (9) << S/4.0*sqrt(2.0*J*S)*sqrt(1-cos(2.0*M_PI*KX))/sqrt(2.0*J*S*(1+cos(2.0*M_PI*KX))+2.0*D*S) << endl;
        
        //cout << endl;
    }
    return 0;
}
