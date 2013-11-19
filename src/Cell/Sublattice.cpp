#include <iostream>
#include "Sublattice.h"

using namespace std;

Sublattice::Sublattice()
{
    this->setName("");
    this->setType("None");
    this->setMoment(1.0,0.0,0.0);
}

void Sublattice::setName(string nameInput)
{
    name = nameInput;
}
    
string Sublattice::getName()
{
    return name;
}

void Sublattice::setType(string typeInput)
{
    type = typeInput;
}

string Sublattice::getType()
{
    return type;
}

void Sublattice::setMoment(double spinInput, double thetaInput, double phiInput)
{
    assert(spinInput > 0.0);
    
    spin = spinInput;
    
    assert(thetaInput >= 0.0 && thetaInput <= M_PI);
    
    while (phiInput > 2.0*M_PI)
    {
        phiInput -= 2.0*M_PI;
    }
    while (phiInput < 0.0)
    {
        phiInput += 2.0*M_PI;
    }
    
    theta = thetaInput;
    phi = phiInput;

    //rotation matrix defined in equation A.1 in J. Phys.: Condens. Matter 21 (2009) 216001
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

double Sublattice::getMoment()
{
    return spin;
}

double Sublattice::getTheta()
{
    return theta;
}

double Sublattice::getPhi()
{
    return phi;
}

Matrix3* Sublattice::getRotationMatrix()
{
    return &rotationMatrix;
}

Matrix3* Sublattice::getInverseMatrix()
{
    return &inverseMatrix;
}

void Sublattice::addAtom(double x, double y, double z)
{
    Vector3 pos(x,y,z);
    position.push_back(pos);
}

Sublattice::Iterator Sublattice::begin()
{
    return Iterator(position.begin());
}

Sublattice::Iterator Sublattice::end()
{
    return Iterator(position.end());
}
