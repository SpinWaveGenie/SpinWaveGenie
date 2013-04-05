//
//  CreateMatrix.h
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __Spin_Wave_Fit__CreateMatrix__
#define __Spin_Wave_Fit__CreateMatrix__

#define _USE_MATH_DEFINES

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
    void CreateMatrix_exchange(double KXP, double KYP, double KZP);
    void CreateMatrix_DMy(double KXP, double KYP, double KZP);
    void CreateMatrix_anis_x();
    void CreateMatrix_anis_z();
    void CreateMatrix_bfield();
    void CreateMatrix_NVO(double KXP, double KYP, double KZP);
    void CreateMatrix_YFeO3(double KXP, double KYP, double KZP);
    void CreateMatrix_AFM(double KXP, double KYP, double KZP);
    void Calc_Eigenvalues();
    void Calc_Weights();
    void Rotation_Matrix();
    //void Rotation_Matrix_NVO();
    void Unique_Solutions();
    void Signif_Solutions(double KXP, double KYP, double KZP);
    std::vector<double> Get_Frequencies();
    std::vector<double> Get_Intensities();

private:
    int M,N;
    //int mod(int K);
    int NU,MI,IM;
    //VectorXd VP;
    std::vector<SW_sublattice> SL;
    std::vector<double> X;
    Eigen::MatrixXcd LN;
    Eigen::MatrixXd CP,CM,CZ;
    Eigen::VectorXd TM,SS,AL;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    std::vector<std::pair<double, std::pair<std::complex<double>,Eigen::VectorXcd> > > eigen;
    Eigen::VectorXd WW;
    Eigen::VectorXd SXX,SYY,SZZ;
    Eigen::VectorXd WP;
    Eigen::VectorXd VP,TXX,TYY,TZZ;
    std::vector<double> VI,SVI;
    Eigen::MatrixXcd XY,XIN;
};

#endif /* defined(__Spin_Wave_Fit__CreateMatrix__) */
