#include "Exch_Interaction.h"
#include "../Cell/AtomIterator.h"
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

Exch_Interaction::Exch_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    this->Update_Interaction(value_in, sl_r_in, sl_s_in, min_in, max_in);
 
}

void Exch_Interaction::Update_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    sl_s = sl_s_in;
    min = min_in;
    max = max_in;
}

vector<string> Exch_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    sl.push_back(sl_s);
    return sl;
}

void Exch_Interaction::calcConstantValues(boost::shared_ptr<Cell> cell)
{
    //find location of r,s
    r= -1;
    s= -1;
    
    M=0;
    for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
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
    
    Sr = cell->getSublattice(sl_r).getMoment();
    Ss = cell->getSublattice(sl_s).getMoment();

    Frs = (*cell->getSublattice(sl_r).getRotationMatrix())*
    (*cell->getSublattice(sl_s).getInverseMatrix());
    
    Fsr = (*cell->getSublattice(sl_s).getRotationMatrix())*
    (*cell->getSublattice(sl_r).getInverseMatrix());
    
    //cout << r << "\t" << s << endl << F << endl;
    //cout << endl;
    
    //cout << "G1= " << G1 << endl;
    //cout << "G2= " << G2 << endl;
    
}

void Exch_Interaction::calcChangingValues(boost::shared_ptr<Cell> cell, Vector3d K)
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
        //cout << K[0]/M_PI << "  " << K[1]/M_PI << "  " << K[2]/M_PI << endl;
        //cout << (*nbr)[0] << "  " << (*nbr)[1] << "  " << (*nbr)[2] << " " << MXI*dot_prod/M_PI << endl;
        //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;

    }
    
    gamma_rs /= z_rs; //force gamma_rs(k=0) = 1.0
    //cout << "z_rs= " << z_rs << endl;
    //cout << "gamma_rs(" << r << "," << s << ")= " << gamma_rs << endl;
}

void Exch_Interaction::checkFirstOrderTerms(boost::shared_ptr<Cell> cell, VectorXcd &elements )
{
        
    complex<double> F1rs(Frs(0,2),Frs(1,2));
    complex<double> F2rs(Frs(2,0),Frs(2,1));
    complex<double> F1sr(Fsr(0,2),Fsr(1,2));
    complex<double> F2sr(Fsr(2,0),Fsr(2,1));
    
    Neighbors neighborList(cell,sl_r,sl_s,min,max);
    AtomIterator nbrBegin = neighborList.begin();
    AtomIterator nbrEnd = neighborList.end();
    z_rs = (double) distance(nbrBegin,nbrEnd);
        
    //elements[r] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(conj(F1rs));
    //elements[s] -= sqrt(Ss)*Sr/(2.0*sqrt(2.0))*z_rs*value*(conj(F2rs));
    //elements[r+M] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(F1rs);
    //elements[s+M] -= sqrt(Ss)*Sr/(2.0*sqrt(2.0))*z_rs*value*(F2rs);
    elements[r] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*conj(F1rs + F2sr);
    elements[r+M] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(F1rs+F2sr);
}

void Exch_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN, int quadrant)
{
    
    complex<double> G1rs = -0.5*complex<double>(Frs(0,0) + Frs(1,1),Frs(1,0)-Frs(0,1));
    
    complex<double> G2rs = -0.5*complex<double>(Frs(0,0) - Frs(1,1),-Frs(1,0)-Frs(0,1));
    
    complex<double> G1sr = -0.5*complex<double>(Fsr(0,0) + Fsr(1,1),Fsr(1,0)-Fsr(0,1));
    
    complex<double> G2sr = -0.5*complex<double>(Fsr(0,0) - Fsr(1,1),-Fsr(1,0)-Fsr(0,1));

    double X = value*sqrt(Sr*Ss);
    switch (quadrant)
    {
        case 0:
            //cout << "G2rs  G2sr" << endl;
            //cout << G2rs << " " << G2sr << endl;
            LN(r,r) += 0.25*z_rs*value*Ss*(Frs(2,2)+Fsr(2,2));
            LN(r,s) += 0.25*z_rs*X*conj(gamma_rs)*(G1rs+conj(G1sr));
            break;
        case 1:
            LN(r,s+M) += 0.25*z_rs*X*conj(gamma_rs)*(conj(G2rs)+conj(G2sr));
            break;
        case 2:
            LN(r+M,s) += 0.25*z_rs*X*conj(gamma_rs)*(G2rs+G2sr);
            break;
        case 3:
            LN(r+M,r+M) += 0.25*z_rs*value*Ss*(Frs(2,2)+Fsr(2,2));
            LN(r+M,s+M) += 0.25*z_rs*X*conj(gamma_rs)*(conj(G1rs)+G1sr);
            break;
        default:
            cout << "error: case must be between 0 and 3" << endl;
            break;
    }
    /*switch (quadrant)
    {
        case 0:
            LN(r,r) += 0.25*z_rs*value*Ss*Frs(2,2);
            LN(r,s) += 0.25*z_rs*X*conj(gamma_rs)*G1rs;
            LN(s,r) += 0.25*z_rs*X*gamma_rs*conj(G1rs);
            LN(s,s) += 0.25*z_rs*value*Sr*Frs(2,2);
            break;
        case 1:
            LN(r,s+M) += 0.25*z_rs*X*conj(gamma_rs)*conj(G2rs);
            LN(s,r+M) += 0.25*z_rs*X*gamma_rs*conj(G2rs);
            break;
        case 2:
            LN(r+M,s) += 0.25*z_rs*X*conj(gamma_rs)*G2rs;
            LN(s+M,r) += 0.25*z_rs*X*gamma_rs*G2rs;
            break;
        case 3:
            LN(r+M,r+M) += 0.25*z_rs*value*Ss*Frs(2,2);
            LN(r+M,s+M) += 0.25*z_rs*X*conj(gamma_rs)*conj(G1rs);
            LN(s+M,r+M) += 0.25*z_rs*X*gamma_rs*G1rs;
            LN(s+M,s+M) += 0.25*z_rs*value*Sr*Frs(2,2);
            break;
        default:
            cout << "error: case must be between 0 and 3" << endl;
            break;
    }*/
}