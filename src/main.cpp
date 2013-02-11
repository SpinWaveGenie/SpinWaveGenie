//
//  main.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#include <iostream>
#include <vector>
#include "SW_Matrix.h"
#include "SW_sublattice.h"

using namespace Eigen;
using namespace std;

int main(int argc, const char * argv[])
{
    double PI = atan(1.0)*4.0;
    vector<double> X(4,0.0);
    vector<SW_sublattice> SL(2);
    double KX;
    SW_Matrix test;
    
    SL[0].set_sublattice(5.0/2.0,PI/2.0 - 0.01098,0.0);
    SL[1].set_sublattice(5.0/2.0,PI/2.0 - 0.01098,PI);
    
    cout << SL[0].get_rot_matrix() << endl;
    cout << SL[1].get_rot_matrix() << endl;
    
    X[0] = 2.112; //2.48;
    X[1] = 0.109;
    X[2] = 4.61e-3;
    X[3] = 1.135e-3;
    
    /*X[0] = -0.1594199E-02;
    X[1] = -0.2192590E+00;
    X[2] = -0.2875218E-01;
    X[3] =  0.3184059E-02;
    X[4] =  0.9540724E-03;
    X[5] =  0.1066680E+01;
    X[6] =  0.3936251E-04;
    X[7] =  0.2611718E+00;
    */
    
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

