//
//  SW_sublattice.h
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 2/6/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __Spin_Wave_Fit__SW_sublattice__
#define __Spin_Wave_Fit__SW_sublattice__

#include <iostream>
#include <Eigen/Dense>
#include <vector>



class SW_sublattice
{
public:
    void set_sublattice(double spin_input, double theta_input, double phi_input);
    void add_neighbors(double type, Eigen::Vector3d position);
    std::vector<double> get_sublattice();
    Eigen::Matrix3d get_rot_matrix();
    Eigen::Matrix3d get_inv_matrix();
private:
    double PI = atan(1.0)*4.0;
    double spin,theta,phi;
    Eigen::Matrix3d rot_matrix,inv_matrix;
};

#endif /* defined(__Spin_Wave_Fit__SW_sublattice__) */
