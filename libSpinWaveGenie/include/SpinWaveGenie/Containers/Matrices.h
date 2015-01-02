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

namespace SpinWaveGenie
{

// Instead of directly importing Vector3d and Matrix3d from the Eigen library, we provide typedefs for them here.
// I've started removing Vector3 from public interfaces because alternatives (such as three doubles) are similarly
// clear, yet don't depend on Eigen.

typedef Eigen::Matrix<double, 3, 1> Vector3;
typedef Eigen::Matrix<double, 3, 3> Matrix3;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXcdRowMajor;
}
#endif /* defined(__spin_wave_genie__Matrices__) */
