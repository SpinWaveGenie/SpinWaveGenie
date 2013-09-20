#include "Sublattice.h"
#include <iostream>
#include "AtomIterator.h"

using namespace std;

void Sublattice::setName(string nameInput)
{
    name = nameInput;
}
    
string Sublattice::getName()
{
    return name;
}

void Sublattice::setMoment(double spinInput, double thetaInput, double phiInput)
{
    spin = spinInput;
    assert(thetaInput <= M_PI && thetaInput >= 0.0);
    theta = thetaInput;
    
    assert(phiInput <= 2.0*M_PI && phiInput >= 0.0);
    phi = phiInput;
    
    angles.push_back(spin);
    angles.push_back(theta);
    angles.push_back(phi);
    
    /*
    rotation matrix defined in equation A.1 in J. Phys.: Condens. Matter 21 (2009) 216001 
    */
    rotationMatrix(0,0) = cos(theta)*cos(phi);
    rotationMatrix(0,1) = cos(theta)*sin(phi);
    rotationMatrix(0,2) = -1.0*sin(theta);
    rotationMatrix(1,0) = -1.0*sin(phi);
    rotationMatrix(1,1) = cos(phi);
    rotationMatrix(1,2) = 0.0;
    rotationMatrix(2,0) = sin(theta)*cos(phi);
    rotationMatrix(2,1) = sin(theta)*sin(phi);
    rotationMatrix(2,2) = cos(theta);
    
    inverseMatrix = rotationMatrix.inverse();
}

vector<double> Sublattice::getMoment()
{
    return angles;
}

Eigen::Matrix3d Sublattice::getRotationMatrix()
{
    return rotationMatrix;
}

Eigen::Matrix3d Sublattice::getInverseMatrix()
{
    return inverseMatrix;
}

void Sublattice::addAtom(double x, double y, double z)
{
    //written to compile with C++03 compiler 
    double tmp[3] = {x,y,z};
    std::vector<double> pos(&tmp[0], &tmp[0]+3);
    position.push_back(pos);
}

AtomIterator Sublattice::begin()
{
    return AtomIterator(position.begin());
}
AtomIterator Sublattice::end()
{
    return AtomIterator(position.end());
}
