//
//  Anisotropy_Interaction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/25/13.
//
//

#include "AnisotropyInteraction.h"
using namespace std;
AnisotropyInteraction::AnisotropyInteraction(double value_in, Vector3 unitVectorIn, string sl_r_in)
{
    this->UpdateInteraction(value_in, unitVectorIn, sl_r_in);
}

Interaction* AnisotropyInteraction::do_clone() const
{
    return new AnisotropyInteraction(*this);
}

void AnisotropyInteraction::UpdateInteraction(double value_in, Vector3 unitVectorIn, string sl_r_in)
{
    value = value_in;
    unitVectorIn.normalize();
    directions = unitVectorIn*unitVectorIn.transpose();
    sl_r = sl_r_in;
}

vector<string> AnisotropyInteraction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    return sl;
}

void AnisotropyInteraction::calcConstantValues(Cell& cell)
{
    complex<double> XI (0.0,1.0);
    //find location of r
    r = cell.getPosition(sl_r);
    M = cell.size();
    double S = cell.getSublattice(sl_r).getMoment();
    Matrix3& inv = cell.getSublattice(sl_r).getInverseMatrix();
    
    //cout << inv << endl;
    
    LNrr = complex<double>(0.0,0.0);
    LNrMr = complex<double>(0.0,0.0);
    LNrrM = complex<double>(0.0,0.0);
    LNrMrM = complex<double>(0.0,0.0);
    
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
            //cout << i << " " << j << " " << directions(i,j) << endl;
            if (directions(i,j) > 1.0e-10)
            {
                double X = value*S*directions(i,j);
                complex<double> beta =    0.5*(inv(i,0)*inv(j,0) - inv(i,1)*inv(j,1) - XI*inv(i,0)*inv(j,1) - XI*inv(i,1)*inv(j,0));
                complex<double> epsilon = 0.5*(inv(i,0)*inv(j,0) + inv(i,1)*inv(j,1) + XI*inv(i,0)*inv(j,1) - XI*inv(i,1)*inv(j,0));
                //cout << "beta= " << beta << " epsilon= " << epsilon << endl;
                LNrr += X*(real(beta) - inv(i,2)*inv(j,2));
                LNrMr += X*epsilon;
                LNrrM += X*conj(epsilon) ;
                LNrMrM += X*(real(beta) - inv(i,2)*inv(j,2));
            }
        }
    //cout << "new implementation" << endl;
    //cout << LNrr << " "<< LNrMr << " " << LNrrM << " " << LNrMrM << endl;
}

void AnisotropyInteraction::checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
    complex<double> XI (0.0,1.0);
    double S = cell.getSublattice(sl_r).getMoment();
    Matrix3& inv = cell.getSublattice(sl_r).getInverseMatrix();
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            //cout << i << " " << j << " " << directions(i,j) << endl;
            if (directions(i,j) > 1.0e-10)
            {
                double X = value*directions(i,j)*sqrt(pow(S,3)/2.0);
                complex<double> zeta = inv(i,2)*inv(j,0) + inv(i,0)*inv(j,2) - XI*(inv(i,2)*inv(j,1) + inv(i,1)*inv(j,2));
                //cout << "zeta= " << zeta << endl;
                elements(r) += X*zeta;
                elements(r+M) += X*conj(zeta);
            }
        }
    }
}

void AnisotropyInteraction::Update_Matrix(Vector3 K, Eigen::MatrixXcd &LN)
{

    LN(r,r) += LNrr;
    LN(r,r+M) += LNrrM;
    LN(r+M,r) += LNrMr;
    LN(r+M,r+M) += LNrMrM;

}
