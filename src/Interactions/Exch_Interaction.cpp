#include "Exch_Interaction.h"
#include <iostream>
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

Exch_Interaction::Exch_Interaction(double value_in, string sl_r_in,string sl_s_in, double min_in, double max_in)
{
    this->Update_Interaction(value_in, sl_r_in, sl_s_in, min_in, max_in);
 
}

Interaction* Exch_Interaction::do_clone() const
{
    return new Exch_Interaction(*this);
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
    vector<string> sl = {sl_r,sl_s};
    return sl;
}

void Exch_Interaction::calcConstantValues(Cell& cell)
{
    r = cell.getPosition(sl_r);
    s = cell.getPosition(sl_s);
    M = cell.size();
    
    Sr = cell.getSublattice(sl_r).getMoment();
    Ss = cell.getSublattice(sl_s).getMoment();

    Frs = cell.getSublattice(sl_r).getRotationMatrix()*
          cell.getSublattice(sl_s).getInverseMatrix();
    
    Fsr = cell.getSublattice(sl_s).getRotationMatrix()*
          cell.getSublattice(sl_r).getInverseMatrix();
    
    neighbors.findNeighbors(cell,sl_r, sl_s, min, max);
    z_rs = neighbors.getNumberNeighbors();
    
    //cout << r << "\t" << s << endl << F << endl;
    //cout << endl;
    
    //cout << "G1= " << G1 << endl;
    //cout << "G2= " << G2 << endl;
    
}

void Exch_Interaction::checkFirstOrderTerms(Cell& cell, VectorXcd &elements )
{
        
    complex<double> F1rs(Frs(0,2),Frs(1,2));
    complex<double> F2rs(Frs(2,0),Frs(2,1));
    complex<double> F1sr(Fsr(0,2),Fsr(1,2));
    complex<double> F2sr(Fsr(2,0),Fsr(2,1));
    
    //elements[r] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(conj(F1rs));
    //elements[s] -= sqrt(Ss)*Sr/(2.0*sqrt(2.0))*z_rs*value*(conj(F2rs));
    //elements[r+M] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(F1rs);
    //elements[s+M] -= sqrt(Ss)*Sr/(2.0*sqrt(2.0))*z_rs*value*(F2rs);
    elements[r] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*conj(F1rs + F2sr);
    elements[r+M] -= sqrt(Sr)*Ss/(2.0*sqrt(2.0))*z_rs*value*(F1rs+F2sr);
}

void Exch_Interaction::Update_Matrix(Vector3d K, MatrixXcd &LN)
{
    
    complex<double> G1rs = -0.5*complex<double>(Frs(0,0) + Frs(1,1),Frs(1,0)-Frs(0,1));
    
    complex<double> G2rs = -0.5*complex<double>(Frs(0,0) - Frs(1,1),-Frs(1,0)-Frs(0,1));
    
    complex<double> G1sr = -0.5*complex<double>(Fsr(0,0) + Fsr(1,1),Fsr(1,0)-Fsr(0,1));
    
    complex<double> G2sr = -0.5*complex<double>(Fsr(0,0) - Fsr(1,1),-Fsr(1,0)-Fsr(0,1));
    

    gamma_rs = neighbors.getGamma(K);
    
    //cout << gamma_rs << endl;
    
    double X = value*sqrt(Sr*Ss);
    

    //cout << "G2rs  G2sr" << endl;
    //cout << G2rs << " " << G2sr << endl;
    LN(r,r) += 0.25*z_rs*value*Ss*(Frs(2,2)+Fsr(2,2));
    LN(r,s) += 0.25*z_rs*X*conj(gamma_rs)*(G1rs+conj(G1sr));
    LN(r,s+M) += 0.25*z_rs*X*conj(gamma_rs)*(conj(G2rs)+conj(G2sr));
    LN(r+M,s) += 0.25*z_rs*X*conj(gamma_rs)*(G2rs+G2sr);
    LN(r+M,r+M) += 0.25*z_rs*value*Ss*(Frs(2,2)+Fsr(2,2));
    LN(r+M,s+M) += 0.25*z_rs*X*conj(gamma_rs)*(conj(G1rs)+G1sr);
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