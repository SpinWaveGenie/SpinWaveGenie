//
//  EigenStuff.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/18/13.
//
//
#ifndef __spin_wave_genie__Matrices__
#define __spin_wave_genie__Matrices__

#include "Eigen/Core"

using Eigen::Dynamic;
using Eigen::RowMajor;

typedef Eigen::Matrix <double, 3, 1> Vector3;
typedef Eigen::Matrix<double, Dynamic, 1> Vector;

typedef Eigen::Matrix <double, 3, 3> Matrix3;
typedef Eigen::Matrix<double, Dynamic, Dynamic> Matrix;
typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> MatrixRowMajor;

#endif /* defined(__spin_wave_genie__Matrices__) */
