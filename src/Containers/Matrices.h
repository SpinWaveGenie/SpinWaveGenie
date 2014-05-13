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

// Vevtor3d and Matrix3d from the eigen library, so we provide typedefs for them here.
// I've started replacing Vector3 from the public interface because alternatives often exist.

typedef Eigen::Matrix <double, 3, 1> Vector3;
typedef Eigen::Matrix <double, 3, 3> Matrix3;

#endif /* defined(__spin_wave_genie__Matrices__) */
