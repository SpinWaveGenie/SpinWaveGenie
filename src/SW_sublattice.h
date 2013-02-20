//
//  SW_sublattice.h
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 2/6/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __Spin_Wave_Fit__SW_sublattice__
#define __Spin_Wave_Fit__SW_sublattice__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>



class SW_sublattice
{
public:
    struct neighbor
    {
        int type;
        Eigen::Vector3d position;
    };
    void set_interactions(std::vector<double> int_input);
    void set_sublattice(double spin_input, double theta_input, double phi_input);
    void add_neighbors(std::vector<int> types, std::vector<Eigen::Vector3d> positions);
    std::vector<neighbor> get_neighbors();
    std::vector<double> get_sublattice();
    Eigen::Matrix3d get_rot_matrix();
    Eigen::Matrix3d get_inv_matrix();
private:
    std::vector<neighbor> neighbors;
    double spin,theta,phi;
    Eigen::Matrix3d rot_matrix,inv_matrix;
    std::vector<double> exchange_interaction;
};

#endif /* defined(__Spin_Wave_Fit__SW_sublattice__) */
