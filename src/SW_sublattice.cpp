//
//  SW_sublattice.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 2/6/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#include "SW_sublattice.h"

using namespace Eigen;
using namespace std;

void SW_sublattice::set_sublattice(double spin_input, double theta_input, double phi_input)
{
    spin = spin_input;
    
    assert(theta_input <= M_PI && theta_input >= 0.0);
    theta = theta_input;
    
    assert(phi_input <= 2.0*M_PI && phi_input >= 0.0);
    phi = phi_input;
    
    rot_matrix(0,0) = cos(theta)*cos(phi);
    rot_matrix(0,1) = cos(theta)*sin(phi);
    rot_matrix(0,2) = -1.0*sin(theta);
    rot_matrix(1,0) = -1.0*sin(phi);
    rot_matrix(1,1) = cos(phi);
    rot_matrix(1,2) = 0.0;
    rot_matrix(2,0) = sin(theta)*cos(phi);
    rot_matrix(2,1) = sin(theta)*sin(phi);
    rot_matrix(2,2) = cos(theta);
    
    inv_matrix = rot_matrix.inverse();
}

void SW_sublattice::add_neighbors(vector<int> input_types, vector<Vector3d> input_positions)
{
    assert(input_types.size() == input_positions.size());
    
    neighbor add;
    for (int i=0;i < input_types.size(); i++)
    {
        add.type = input_types[i];
        add.position = input_positions[i];
        neighbors.push_back(add);
    }
    
}

vector<double> SW_sublattice::get_sublattice()
{
    vector<double> angles;
    angles.push_back(spin);
    angles.push_back(theta);
    angles.push_back(phi);
    
    //neighbor_list =
    
    return angles;
}

vector<Vector3d> SW_sublattice::get_neighbors(int type)
{
    vector<Vector3d> requested_neighbors;
    
    for(int i=0; i<neighbors.size();i++)
    {
        if (neighbors[i].type == type)
            requested_neighbors.push_back(neighbors[i].position);
    }
    
    return requested_neighbors;
}

Eigen::Matrix3d SW_sublattice::get_rot_matrix()
{
    return rot_matrix;
}

Eigen::Matrix3d SW_sublattice::get_inv_matrix()
{
    return inv_matrix;
}
