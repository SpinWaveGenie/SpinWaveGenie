#include "ExchangeInteractionSameSublattice.h"
#include <iostream>
#include "../Cell/Neighbors.h"

using namespace std;
using namespace Eigen;

namespace SpinWaveGenie
{

ExchangeInteractionSameSublattice::ExchangeInteractionSameSublattice(string name_in, double value_in, string sl_r_in, double min_in, double max_in)
{
    name = name_in;
    this->updateInteraction(value_in, sl_r_in, min_in, max_in);
}

Interaction* ExchangeInteractionSameSublattice::do_clone() const
{
    return new ExchangeInteractionSameSublattice(*this);
}

void ExchangeInteractionSameSublattice::updateInteraction(double value_in, string sl_r_in, double min_in, double max_in)
{
    value = value_in;
    sl_r = sl_r_in;
    min = min_in;
    max = max_in;
}

const string& ExchangeInteractionSameSublattice::getName()
{
    return name;
}

void ExchangeInteractionSameSublattice::updateValue(double value_in)
{
    value = value_in;
}

vector<string> ExchangeInteractionSameSublattice::sublattices() const
{
    vector<string> sl = {sl_r,sl_r};
    return sl;
}

void ExchangeInteractionSameSublattice::calcConstantValues(Cell& cell)
{
    //cout << "cell check(calcConstantValues): " << cell.begin()->getName() << endl;
    r = cell.getPosition(sl_r);
    M = cell.size();
    
    //cout << "r: " << r << " " << ", M: " << M << endl;
    double Sr = cell.getSublattice(sl_r).getMoment();

    //cout << "Sr: " << Sr << endl;

    Matrix3 Frs = cell.getSublattice(sl_r).getRotationMatrix()*
          cell.getSublattice(sl_r).getInverseMatrix();
       
    neighbors.findNeighbors(cell,sl_r, sl_r, min, max);
    double z_rs = neighbors.size();
    //cout << "z_rs = " << z_rs << endl;

    //cout << "cell check(calcConstantValues): " << cell2.begin()->getName() << endl;

    complex<double> G1rs = -0.5*complex<double>(Frs(0,0) + Frs(1,1),Frs(1,0)-Frs(0,1));
    complex<double> G2rs = -0.5*complex<double>(Frs(0,0) - Frs(1,1),-Frs(1,0)-Frs(0,1));
    
    LNrr = 0.25*z_rs*value*Sr*(Frs(2,2)+Frs(2,2));
    LNrs = 0.25*z_rs*value*Sr*(G1rs+conj(G1rs));
    LNrsM = 0.25*z_rs*value*Sr*conj(G2rs+G2rs);
    
    //cout << LNrr << " " << LNrs << " " << LNrsM << endl;
}

void ExchangeInteractionSameSublattice::calculateEnergy(Cell& cell, double &energy)
{
    r = cell.getPosition(sl_r);
    neighbors.findNeighbors(cell,sl_r, sl_r, min, max);
    double z_rs = neighbors.size();
    double Sr = cell[r].getMoment();
    
    energy -= value*z_rs*Sr*Sr;
}

void ExchangeInteractionSameSublattice::calculateFirstOrderTerms(Cell& cell, VectorXcd &elements )
{
    //first order terms are 0.0
}

void ExchangeInteractionSameSublattice::updateMatrix(Vector3d K, MatrixXcd &LN)
{
    gamma_rs = neighbors.getGamma(K);
    //cout << value << " " << sl_r << " " << sl_r << " " << gamma_rs << endl;

    //cout << "number of neighbors: " << neighbors.size() << endl;
    //cout << gamma_rs <<  " " << conj(gamma_rs) << endl;
    //cout << LNrr/value << endl;
    
    LN(r,r) += LNrr + LNrs*conj(gamma_rs);
    LN(r,r+M) += LNrsM*conj(gamma_rs);
    LN(r+M,r) += conj(LNrsM)*conj(gamma_rs);
    LN(r+M,r+M) += LNrr + conj(LNrs)*conj(gamma_rs);
    
}

}