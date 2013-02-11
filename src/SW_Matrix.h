//
//  CreateMatrix.h
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __Spin_Wave_Fit__CreateMatrix__
#define __Spin_Wave_Fit__CreateMatrix__

#include <iostream>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <ctime>
#include "SW_sublattice.h"

class SW_Matrix
{
public:
    void set_parr(std::vector<SW_sublattice> SLin, std::vector<double>& Xin);
    std::vector<double> get_parr();
    void CreateMatrix_NVO(double KXP, double KYP, double KZP);
    void CreateMatrix_YFeO3(double KXP, double KYP, double KZP);
    void CreateMatrix_AFM(double KXP, double KYP, double KZP);
    void Calc_Eigenvalues();
    void Calc_Weights();
    void Rotation_Matrix();
    //void Rotation_Matrix_NVO();
    void Unique_Solutions();
    void Signif_Solutions(double KXP, double KYP, double KZP);
private:
    /*const int P = 10;
    const int M = 20;
    const int N = 40;*/
    /*const int P = 1;
    const int M = 2;
    const int N = 4;*/
    int M,N;
    //int mod(int K);
    int NU,MI,IM;
    //VectorXd VP;
    std::vector<SW_sublattice> SL;
    std::vector<double> X;
    double PI = atan(1.0)*4.0;
    Eigen::MatrixXcd LN;
    Eigen::MatrixXd CP,CM,CZ;
    Eigen::VectorXd TM,SS,AL;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    std::vector<std::pair<double, std::pair<std::complex<double>,Eigen::VectorXcd> > > eigen;
    Eigen::VectorXd WW;
    Eigen::VectorXd SXX,SYY,SZZ;
    Eigen::VectorXd WP;
    Eigen::VectorXd VP,TXX,TYY,TZZ;
    Eigen::VectorXd VI,SVI;
    Eigen::MatrixXcd XY,XIN;
};

#endif /* defined(__Spin_Wave_Fit__CreateMatrix__) */
