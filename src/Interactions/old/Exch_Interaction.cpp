#include "Exch_Interaction.h"

using namespace std;
using namespace Eigen;

Exch_Interaction::Exch_Interaction()
{
    
}

void Exch_Interaction::Add_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    Exch_Parameters in;
    in.value = value_in;
    in.sl_r = sl_r_in;
    in.sl_s = sl_s_in;
    in.min = min_in;
    in.max = max_in;
    exch_array.push_back(in);
}

void Exch_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN)
{
    //K[0] = KXP*2.0*M_PI;
    //K[1] = KYP*2.0*M_PI;
    //K[2] = KZP*2.0*M_PI; 
    
    std::vector<Exch_Parameters>::iterator iter;
    for(iter=exch_array.begin();iter!=exch_array.end();iter++)
    {
        //find location of r,s
        int r= -1;
        int s= -1;
    
        CellIter sl(cell);
        int M=0;
        for(sl.First();!sl.IsDone();sl.Next())
        {
            if ( iter->sl_r == sl.CurrentItem()->get_name())
                r = M;
            if ( iter->sl_s == sl.CurrentItem()->get_name())
                s = M;
            M++;
        }
        assert(r!=-1 && s!=-1);
           
        double S = cell->get_sublattice(iter->sl_r)->get_moment()[0];
        
        Matrix3d F;
        F = cell->get_sublattice(iter->sl_r)->get_rot_matrix()*
            cell->get_sublattice(iter->sl_s)->get_inv_matrix();
    
        //cout << r << "\t" << s << endl << F << endl;
        //cout << endl;

        boost::shared_ptr<Sublattice> sl_s = cell->get_sublattice(iter->sl_r);
        boost::shared_ptr<Sublattice> sl_r = cell->get_sublattice(iter->sl_s);
        NeighborIter nbr(cell,sl_s,sl_r,iter->min,iter->max);

        double z_rs;
        z_rs = 0;
        complex<double> gamma_rs (0.0,0.0);
        for(nbr.First();!nbr.IsDone();nbr.Next())
        {
            double dot_prod = K.dot(nbr.CurrentItem());
            //cout << nbr.CurrentItem().transpose() << endl;
            gamma_rs += complex<double> (cos(dot_prod),-1.0*sin(dot_prod));
            z_rs = z_rs + 1.0;
        }
        
        gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
        //cout << "z_rs= " << z_rs << endl;
        //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;
    
        complex<double> G1 (F(0,0) + F(1,1),F(1,0)-F(0,1));
        G1 *= -0.5;
    
        complex<double> G2 (F(0,0) - F(1,1),-F(1,0)-F(0,1));
        G2 *= -0.5;
        
        double X = iter->value;
                
        //cout << "G1= " << G1 << endl;
	    //cout << "G2= " << G2 << endl;

        LN(r,r) += 0.25*z_rs*X*S*F(2,2);
        LN(s,s) += 0.25*z_rs*X*S*F(2,2);
        LN(r+M,r+M) += 0.25*z_rs*X*S*F(2,2);
        LN(s+M,s+M) += 0.25*z_rs*X*S*F(2,2);
    
        LN(r,s) += 0.25*z_rs*X*S*gamma_rs*G1;
        LN(r+M,s+M) += 0.25*z_rs*X*S*conj(gamma_rs)*conj(G1);
        LN(s,r) += 0.25*z_rs*X*S*gamma_rs*conj(G1);
        LN(s+M,r+M) += 0.25*z_rs*X*S*conj(gamma_rs)*G1;
                    
        LN(r,s+M) += 0.25*z_rs*X*S*conj(gamma_rs)*conj(G2);
        LN(s,r+M) += 0.25*z_rs*X*S*gamma_rs*conj(G2);
                    
        LN(r+M,s) += 0.25*z_rs*X*S*conj(gamma_rs)*G2;
        LN(s+M,r) += 0.25*z_rs*X*S*gamma_rs*G2;
    }
}