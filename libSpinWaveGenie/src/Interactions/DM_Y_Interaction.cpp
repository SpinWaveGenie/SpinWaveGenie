#include "SpinWaveGenie/Interactions/DM_Y_Interaction.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

using namespace std;
using namespace Eigen;

namespace SpinWaveGenie
{

DM_Y_Interaction::DM_Y_Interaction(string name_in, double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    name = name_in;
    this->updateInteraction(value_in, sl_r_in, sl_s_in, min_in, max_in);
}

Interaction* DM_Y_Interaction::do_clone() const
{
    return new DM_Y_Interaction(*this);

}

void DM_Y_Interaction::updateInteraction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    sl_s = sl_s_in;
    min = min_in;
    max = max_in;
}

const string& DM_Y_Interaction::getName()
{
    return name;
}

void DM_Y_Interaction::updateValue(double value_in)
{
    value = value_in;
}

vector<string> DM_Y_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    sl.push_back(sl_s);
    return sl;
}

void DM_Y_Interaction::calcConstantValues(Cell& cell)
{
    r = cell.getPosition(sl_r);
    s = cell.getPosition(sl_s);
    M = cell.size();
    assert(r!=-1 && s!=-1);
    
    //cout << r << "\t" << s << endl << F << endl;
    //cout << endl;
    
    double X = value;
    double S = cell.getSublattice(sl_r).getMoment();
    double theta_r = cell.getSublattice(sl_r).getTheta();
    double phi_r = cell.getSublattice(sl_r).getPhi();
    double theta_s = cell.getSublattice(sl_s).getTheta();
    double phi_s = cell.getSublattice(sl_s).getPhi();
    
    value0 = -0.5*X*S*(sin(theta_r)*cos(theta_s)*cos(phi_r) - cos(theta_r)*sin(theta_s)*cos(phi_s));
    value1 = -0.25*X*S*(cos(theta_r)*sin(theta_s)*cos(phi_r)-sin(theta_r)*cos(theta_s)*cos(phi_s));
    value2 = -0.25*X*S*sin(theta_r)*sin(phi_s);
    value3 = -0.25*X*S*sin(theta_s)*sin(phi_r);
    
    neighbors.findNeighbors(cell,sl_r, sl_s, min, max);
    z_rs = neighbors.size();
    //cout << value0 << " " << value1 << " " << value2 << " " << value3 << " " << endl;
}

void DM_Y_Interaction::calculateEnergy(Cell& cell, double &energy)
{
    
}

void DM_Y_Interaction::calculateFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
}


void DM_Y_Interaction::updateMatrix(Vector3d K, MatrixXcd &LN)
{
    complex<double> XI (0.0,1.0);
    gamma_rs = neighbors.getGamma(K);
    //cout << value0 << " " << value1 << " " << value2 << " " << value3 << " " << endl;
    //cout << z_rs << " " << endl;
    LN(r,r) += z_rs*value0;
    LN(r,s) += z_rs*conj(gamma_rs)*(value1 - XI*value2 - XI*value3);
    LN(s,r) += z_rs*gamma_rs*(value1 + XI*value2 + XI*value3);
    LN(s,s) += z_rs*value0;
    LN(r,s+M) += z_rs*conj(gamma_rs)*(value1 + XI*value2 - XI*value3);
    LN(s,r+M) += z_rs*gamma_rs*(value1 + XI*value2 - XI*value3);
    LN(r+M,s) += z_rs*conj(gamma_rs)*(value1 - XI*value2 + XI*value3);
    LN(s+M,r) += z_rs*gamma_rs*(value1 - XI*value2 + XI*value3);
    LN(r+M,r+M) += z_rs*value0;
    LN(r+M,s+M) += z_rs*conj(gamma_rs)*(value1 + XI*value2 + XI*value3);
    LN(s+M,s+M) += z_rs*value0;
    LN(s+M,r+M) += z_rs*gamma_rs*(value1 - XI*value2 - XI*value3);

}

}