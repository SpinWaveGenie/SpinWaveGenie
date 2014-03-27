#include "ExchangeInteractionSameSublattice.h"
#include <iostream>
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

ExchangeInteractionSameSublattice::ExchangeInteractionSameSublattice(string name_in, double value_in, string sl_r_in, double min_in, double max_in)
{
    name = name_in;
    this->Update_Interaction(value_in, sl_r_in, min_in, max_in);
}

Interaction* ExchangeInteractionSameSublattice::do_clone() const
{
    return new ExchangeInteractionSameSublattice(*this);
}

void ExchangeInteractionSameSublattice::Update_Interaction(double value_in, string sl_r_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    min = min_in;
    max = max_in;
}

string ExchangeInteractionSameSublattice::getName()
{
    return name;
}

void ExchangeInteractionSameSublattice::updateValue(double value_in)
{
    value = value_in;
}

vector<string> ExchangeInteractionSameSublattice::sublattices() const
{
    vector<string> sl = {sl_r,sl_s};
    return sl;
}

void ExchangeInteractionSameSublattice::calcConstantValues(Cell& cell)
{
    r = cell.getPosition(sl_r);
    M = cell.size();
    
    double Sr = cell.getSublattice(sl_r).getMoment();

    Matrix3 Frs = cell.getSublattice(sl_r).getRotationMatrix()*
          cell.getSublattice(sl_s).getInverseMatrix();
    
    neighbors.findNeighbors(cell,sl_r, sl_s, min, max);
    double z_rs = neighbors.getNumberNeighbors();
    
    complex<double> F1rs(Frs(0,2),Frs(1,2));
    complex<double> F2sr(Frs(2,0),Frs(2,1));
    
    LNr = -1.0*sqrt(Sr)*Sr/(2.0*sqrt(2.0))*z_rs*value*conj(F1rs + F2sr);

    complex<double> G1rs = -0.5*complex<double>(Frs(0,0) + Frs(1,1),Frs(1,0)-Frs(0,1));
    complex<double> G2rs = -0.5*complex<double>(Frs(0,0) - Frs(1,1),-Frs(1,0)-Frs(0,1));
    
    LNrr = 0.25*z_rs*value*Sr*(Frs(2,2)+Frs(2,2));
    LNrs = 0.25*z_rs*value*Sr*(G1rs+conj(G1rs));
    LNrsM = 0.25*z_rs*value*Sr*conj(G2rs+G2rs);
}

void ExchangeInteractionSameSublattice::checkFirstOrderTerms(Cell& cell, VectorXcd &elements )
{
    elements[r] += LNr;
    elements[r+M] += conj(LNr);
}

void ExchangeInteractionSameSublattice::Update_Matrix(Vector3d K, MatrixXcd &LN)
{
    gamma_rs = neighbors.getGamma(K);
    
    LN(r,r) += LNrr + LNrs*conj(gamma_rs);
    LN(r,r+M) += LNrsM*conj(gamma_rs);
    LN(r+M,r) += conj(LNrsM)*conj(gamma_rs);
    LN(r+M,r+M) += LNrr + conj(LNrs)*conj(gamma_rs);
    
}