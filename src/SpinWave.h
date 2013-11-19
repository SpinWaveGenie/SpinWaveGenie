//
//  SpinWave.h
//  Spin Wave Genie
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __SpinWave_H__
#define __SpinWave_H__

#define _USE_MATH_DEFINES
#include <iostream>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <ctime>
#include "Cell.h"
#include "Cell/Neighbors.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include "Interaction.h"

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixXcdRowMajor;

struct results
{
    double weight;
    long index;
    // < operator reversed for sort!!!!
    bool operator<( const results& val ) const {
    	return weight > val.weight;
        }
};

//! SpinWave Class
/*!
The SpinWave Class stores the "L" matrix and calculates
the spin wave frequencies and intensities.
*/
class SpinWave
{
public:
    //! Use SW_Builder to generate SpinWave instance
    friend class SW_Builder;
    SpinWave();
    SpinWave(Cell& cell_in, boost::ptr_vector<Interaction> interactions_in);
    //!
    void Set_Kpoint(double KX, double KY, double KZ);
    void Clear_Matrix();
    Eigen::VectorXcd checkFirstOrderTerms();
    void createMatrix(double KX,double KY,double KZ);
    void Calc();
    std::vector<double> Get_Frequencies();
    std::vector<double> Get_Intensities();
    //Eigen::VectorXd SpinWave::Get_Evector(double E_min, double E_max, double E_points);
    //void save(std::string filename);
    //Eigen::MatrixXcd get_Matrix();
private:
    double KXP,KYP,KZP;
    Cell cell;
    Eigen::Vector3d kpoint;
    void Calc_Eigenvalues();
    void Calc_Weights();
    void Calc_Intensities();
    void Unique_Solutions();
    void Signif_Solutions();
    size_t M,N;
    int NU,MI,IM;
    Eigen::MatrixXcd LN;
    Eigen::VectorXd SS;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    Eigen::VectorXd WW; // want to get rid of this
    std::vector<double> VI,SVI; 
    Eigen::MatrixXcd XY,XIN;
    boost::ptr_vector<Interaction> interactions;

};

#endif /* defined(__SpinWave__) */
