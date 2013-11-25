#include "Anis_Z_Interaction.h"
#include <iostream>
using namespace std;
using namespace Eigen;


Anis_Z_Interaction::Anis_Z_Interaction(double value_in, string sl_r_in)
{
    this->Update_Interaction(value_in, sl_r_in);
}

Interaction* Anis_Z_Interaction::do_clone() const
{
    return new Anis_Z_Interaction(*this);
}

void Anis_Z_Interaction::Update_Interaction(double value_in, string sl_r_in)
{
    value = value_in;
    sl_r = sl_r_in;
}

vector<string> Anis_Z_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    return sl;
}

void Anis_Z_Interaction::calcConstantValues(Cell& cell)
{
    r = cell.getPosition(sl_r);
    M = cell.size();

    double S = cell.getSublattice(sl_r).getMoment();
    double theta = cell.getSublattice(sl_r).getTheta();
    double X = value;
    
    LNrr = -0.5*X*S*(1.0-3.0*pow(cos(theta),2));
    LNrrM = -0.5*X*S*pow(sin(theta),2);
    LNrMr = -0.5*X*S*pow(sin(theta),2);
    LNrMrM = -0.5*X*S*(1.0-3.0*pow(cos(theta),2));
    
    cout << "working Z implementation" << endl;
    cout << LNrr << " "<< LNrMr << " " << LNrrM << " " << LNrMrM << endl;
}

void Anis_Z_Interaction::checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
    double S = cell.getSublattice(sl_r).getMoment();
    double theta = cell.getSublattice(sl_r).getTheta();
    
    elements[r] += sqrt(pow(S,3)/2.0)*value*sin(2.0*theta);
    elements[r+M] += sqrt(pow(S,3)/2.0)*value*sin(2.0*theta);
}

void Anis_Z_Interaction::Update_Matrix(Vector3d K, MatrixXcd &LN)
{
    LN(r,r) += LNrr;
    LN(r+M,r) += LNrMr;
    LN(r,r+M) += LNrrM;
    LN(r+M,r+M) += LNrMrM;
    
}