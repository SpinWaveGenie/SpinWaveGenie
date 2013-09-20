#include "DM_Y_Interaction.h"
#include "../Cell/AtomIterator.h"
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

void DM_Y_Interaction::Add_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    sl_s = sl_s_in;
    min = min_in;
    max = max_in;
}

void DM_Y_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN)
{
    double tmp;
    complex<double> tmp1,tmp2,tmp3;
    complex<double> XI (0.0,1.0);

    //find location of r,s
    int r= -1;
    int s= -1;
    
    int M=0;
    for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
    {
        //cout << sl.CurrentItem()->get_name() << endl;
        //cout << iter->sl_r << endl;
        //cout << iter->sl_s << endl;
        if ( sl_r == (*sl)->getName())
            r = M;
        if ( sl_s == (*sl)->getName())
            s = M;
        M++;
    }
    assert(r!=-1 && s!=-1);
                      
    //cout << r << "\t" << s << endl << F << endl;
    //cout << endl;

    boost::shared_ptr<Sublattice> sl_sp = cell->get_sublattice(sl_r);
    boost::shared_ptr<Sublattice> sl_rp = cell->get_sublattice(sl_s);
    
    Neighbors neighbor_list(cell,sl_rp,sl_sp,min,max);
    
    AtomIterator nbrBegin = neighbor_list.begin();
    AtomIterator nbrEnd = neighbor_list.end();
    double z_rs = (double) distance(nbrBegin,nbrEnd);
    
    complex<double> MXI (0.0,-1.0);
    complex<double> gamma_rs (0.0,0.0);
    for(AtomIterator nbr=nbrBegin;nbr!=nbrEnd;++nbr)
    {
        double dot_prod = K[0]*(*nbr)[0] + K[1]*(*nbr)[1] + K[2]*(*nbr)[2];
        gamma_rs += exp(MXI*dot_prod);
    }
        
    gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
    //cout << "z_rs= " << z_rs << endl;
    //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;
            
    double X = value;
    double S = cell->get_sublattice(sl_r)->getMoment()[0];
    double theta_r = cell->get_sublattice(sl_r)->getMoment()[1];
    double phi_r = cell->get_sublattice(sl_r)->getMoment()[2];
    double theta_s = cell->get_sublattice(sl_s)->getMoment()[1];
    double phi_s = cell->get_sublattice(sl_s)->getMoment()[2];
    
    tmp = 0.5*X*S*z_rs*(sin(theta_r)*cos(theta_s)*cos(phi_r) - cos(theta_r)*sin(theta_s)*cos(phi_s));
    LN(r,r) -= tmp;
    LN(r+M,r+M) -= tmp;
    LN(s,s) -= tmp;
    LN(s+M,s+M) -= tmp;
        
    tmp1 = cos(theta_r)*sin(theta_s)*cos(phi_r)-sin(theta_r)*cos(theta_s)*cos(phi_s);
    tmp2 = sin(theta_r)*sin(phi_s);
    tmp3 = sin(theta_s)*sin(phi_r);
        
    //cout << "r= " << r << ", s= " << s << endl;
    //cout << tmp1 << '\t' << tmp2 << '\t' << tmp3 << endl;
    
    LN(r,s) -= 0.25*X*S*z_rs*conj(gamma_rs)*(tmp1 - XI*tmp2 - XI*tmp3);
    LN(s,r) -= 0.25*X*S*z_rs*gamma_rs*(tmp1 + XI*tmp2 + XI*tmp3);
    LN(r+M,s+M) -= 0.25*X*S*z_rs*conj(gamma_rs)*(tmp1 + XI*tmp2 + XI*tmp3);
    LN(s+M,r+M) -= 0.25*X*S*z_rs*gamma_rs*(tmp1 - XI*tmp2 - XI*tmp3);
    
    LN(r+M,s) -= 0.25*X*S*z_rs*conj(gamma_rs)*(tmp1 - XI*tmp2 + XI*tmp3);
    LN(s+M,r) -= 0.25*X*S*z_rs*gamma_rs*(tmp1 - XI*tmp2 + XI*tmp3);
        
    LN(r,s+M) -= 0.25*X*S*z_rs*conj(gamma_rs)*(tmp1 + XI*tmp2 - XI*tmp3);
    LN(s,r+M) -= 0.25*X*S*z_rs*gamma_rs*(tmp1 + XI*tmp2 - XI*tmp3);
}