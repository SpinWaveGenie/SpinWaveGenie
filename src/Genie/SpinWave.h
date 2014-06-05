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
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include "Containers/Results.h"
#include "Interactions/Interaction.h"
#include "Genie/MagneticFormFactor.h"
#include "Containers/Results.h"


typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixXcdRowMajor;

struct results
{
    double weight;
    long index;
};

//! SpinWave Class
/*!
The SpinWave Class stores the "L" matrix and calculates
the spin wave frequencies and intensities.
*/
class SpinWave
{
public:
    SpinWave();
    //! Use SpinWaveBuilder to generate SpinWave instance
    friend class SpinWaveBuilder;
    SpinWave(Cell& cell_in, boost::ptr_vector<Interaction> interactions_in);
    //!
    void setKPoint(double KX, double KY, double KZ);
    void clearMatrix();
    const Cell& getCell() const;
    void createMatrix(double KX,double KY,double KZ);
    void calculate();
    Results getPoints();
private:
    double KXP,KYP,KZP;
    Cell cell;
    void calculateEigenvalues();
    void calculateWeights();
    void calculateIntensities();
    size_t M,N;
    int NU,MI,IM;
    Eigen::MatrixXcd LN;
    Eigen::VectorXd SS;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    Eigen::VectorXd WW; // want to get rid of this
    Results VI;
    Eigen::MatrixXcd XY,XIN;
    boost::ptr_vector<Interaction> interactions;
    MagneticFormFactor formFactor;
};

#endif /* defined(__SpinWave__) */
