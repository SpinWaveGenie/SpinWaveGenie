#include "DM_Y_Interaction.h"
#include "../Cell/AtomIterator.h"
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

DM_Y_Interaction::DM_Y_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    this->Update_Interaction(value_in, sl_r_in, sl_s_in, min_in, max_in);
}

Interaction* DM_Y_Interaction::do_clone() const
{
    return new DM_Y_Interaction(*this);

}

void DM_Y_Interaction::Update_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    sl_s = sl_s_in;
    min = min_in;
    max = max_in;
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
    //find location of r,s
    r= -1;
    s= -1;
    M=0;
    for (SublatticeIterator sl=cell.begin(); sl!=cell.end(); ++sl)
    {
        //cout << sl.CurrentItem()->get_name() << endl;
        //cout << iter->sl_r << endl;
        //cout << iter->sl_s << endl;
        if ( sl_r == sl->getName())
            r = M;
        if ( sl_s == sl->getName())
            s = M;
        M++;
    }
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
    
    //cout << value0 << " " << value1 << " " << value2 << " " << value3 << " " << endl;
    
}

void DM_Y_Interaction::calcChangingValues(Cell& cell, Vector3d K)
{
    Neighbors neighborList(cell,sl_r,sl_s,min,max);
    
    AtomIterator nbrBegin = neighborList.begin();
    AtomIterator nbrEnd = neighborList.end();
    z_rs = (double) distance(nbrBegin,nbrEnd);
    
    complex<double> MXI (0.0,-1.0);
    gamma_rs = complex<double>(0.0,0.0);
    for(AtomIterator nbr=nbrBegin;nbr!=nbrEnd;++nbr)
    {
        double dot_prod = K[0]*(*nbr)[0] + K[1]*(*nbr)[1] + K[2]*(*nbr)[2];
        gamma_rs += exp(MXI*dot_prod);
    }
    
    gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
    //cout << "z_rs= " << z_rs << endl;
    //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;
    
}

void DM_Y_Interaction::checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements)
{
}


void DM_Y_Interaction::Update_Matrix(Vector3d K,Cell& cell, MatrixXcd &LN, int quadrant)
{
    complex<double> XI (0.0,1.0);
    //cout << value0 << " " << value1 << " " << value2 << " " << value3 << " " << endl;
    //cout << z_rs << " " << endl;
    switch (quadrant)
    {
        case 0:
            LN(r,r) += z_rs*value0;
            LN(r,s) += z_rs*conj(gamma_rs)*(value1 - XI*value2 - XI*value3);
            LN(s,r) += z_rs*gamma_rs*(value1 + XI*value2 + XI*value3);
            LN(s,s) += z_rs*value0;
            break;
        case 1:
            LN(r,s+M) += z_rs*conj(gamma_rs)*(value1 + XI*value2 - XI*value3);
            LN(s,r+M) += z_rs*gamma_rs*(value1 + XI*value2 - XI*value3);
            break;
        case 2:
            LN(r+M,s) += z_rs*conj(gamma_rs)*(value1 - XI*value2 + XI*value3);
            LN(s+M,r) += z_rs*gamma_rs*(value1 - XI*value2 + XI*value3);
            break;
        case 3:
            LN(r+M,r+M) += z_rs*value0;
            LN(r+M,s+M) += z_rs*conj(gamma_rs)*(value1 + XI*value2 + XI*value3);
            LN(s+M,s+M) += z_rs*value0;
            LN(s+M,r+M) += z_rs*gamma_rs*(value1 - XI*value2 - XI*value3);
            break;
        default:
            //cout << "error: case must be between 0 and 3" << endl;
            break;
    }
}