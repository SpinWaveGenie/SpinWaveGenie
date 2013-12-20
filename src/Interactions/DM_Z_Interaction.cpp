#include "DM_Z_Interaction.h"
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

DM_Z_Interaction::DM_Z_Interaction(string name_in, double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    name = name_in;
    this->Update_Interaction(value_in, sl_r_in, sl_s_in, min_in, max_in);
}

Interaction* DM_Z_Interaction::do_clone() const
{
    return new DM_Z_Interaction(*this);

}

void DM_Z_Interaction::Update_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    sl_s = sl_s_in;
    min = min_in;
    max = max_in;
}

vector<string> DM_Z_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    sl.push_back(sl_s);
    return sl;
}

string DM_Z_Interaction::getName()
{
    return name;
}

void DM_Z_Interaction::updateValue(double value_in)
{
    value = value_in;
}

void DM_Z_Interaction::calcConstantValues(Cell& cell)
{
    r = cell.getPosition(sl_r);
    s = cell.getPosition(sl_s);
    M = cell.size();
    
    //cout << r << "\t" << s << endl << F << endl;
    //cout << endl;
    
    double X = value;
    double S = cell.getSublattice(sl_r).getMoment();
    double theta_r = cell.getSublattice(sl_r).getTheta();
    double phi_r = cell.getSublattice(sl_r).getPhi();
    double theta_s = cell.getSublattice(sl_s).getTheta();
    double phi_s = cell.getSublattice(sl_s).getPhi();
    
    tmp0 = 0.5*X*S*sin(theta_r)*sin(theta_s)*sin(phi_r-phi_s);
    tmp1 = 0.25*X*S*cos(theta_r)*cos(theta_s)*sin(phi_r-phi_s);
    tmp2 = 0.25*X*S*cos(theta_r)*cos(phi_r-phi_s);
    tmp3 = 0.25*X*S*cos(theta_s)*cos(phi_r-phi_s);
    tmp4 = 0.25*X*S*sin(phi_r-phi_s);
    
    neighbors.findNeighbors(cell,sl_r, sl_s, min, max);
    z_rs = neighbors.getNumberNeighbors();
}


void DM_Z_Interaction::checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
}


void DM_Z_Interaction::Update_Matrix(Vector3d K, MatrixXcd &LN)
{
    complex<double> XI (0.0,1.0);
    gamma_rs = neighbors.getGamma(K);

    LN(r,r) -= z_rs*tmp0;
    LN(r,s) -= z_rs*conj(gamma_rs)*(tmp1 - XI*tmp2 - XI*tmp3 - tmp4);
    LN(s,r) -= z_rs*gamma_rs*(-tmp1 + XI*tmp2 + XI*tmp3 - tmp4);
    LN(s,s) -= z_rs*tmp0;
    LN(r,s+M) -= z_rs*conj(gamma_rs)*(-tmp1 + XI*tmp2 - XI*tmp3 + tmp4);
    LN(s,r+M) -= z_rs*gamma_rs*(-tmp1 + XI*tmp2 - XI*tmp3 + tmp4);
    LN(r+M,s) -= z_rs*conj(gamma_rs)*(-tmp1 - XI*tmp2 + XI*tmp3 + tmp4);
    LN(s+M,r) -= z_rs*gamma_rs*(-tmp1 - XI*tmp2 + XI*tmp3 + tmp4);
    LN(r+M,r+M) -= z_rs*tmp0;
    LN(r+M,s+M) -= z_rs*conj(gamma_rs)*(-tmp1 + XI*tmp2 + XI*tmp3 - tmp4);
    LN(s+M,s+M) -= z_rs*tmp0;
    LN(s+M,r+M) -= z_rs*gamma_rs*(-tmp1 - XI*tmp2 - XI*tmp3 - tmp4);

}