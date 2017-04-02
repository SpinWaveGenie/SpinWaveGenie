#include "SpinWaveGenie/Interactions/ExchangeInteractionSameSublattice.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Memory.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace Eigen;

namespace SpinWaveGenie
{

ExchangeInteractionSameSublattice::ExchangeInteractionSameSublattice(const string &name_in, double value_in,
                                                                     const string &sl_r_in, double min_in,
                                                                     double max_in)
    : name(name_in), r(0), s(0), M(0)
{
  this->updateInteraction(value_in, sl_r_in, min_in, max_in);
}

std::unique_ptr<Interaction> ExchangeInteractionSameSublattice::clone() const
{
  return memory::make_unique<ExchangeInteractionSameSublattice>(*this);
}

void ExchangeInteractionSameSublattice::updateInteraction(double value_in, const string &sl_r_in, double min_in,
                                                          double max_in)
{
  value = value_in;
  sl_r = sl_r_in;
  min = min_in;
  max = max_in;
}

const string &ExchangeInteractionSameSublattice::getName() const { return name; }

void ExchangeInteractionSameSublattice::updateValue(double value_in) { value = value_in; }

std::array<std::string, 2> ExchangeInteractionSameSublattice::sublattices() const { return {{sl_r, sl_r}}; }

void ExchangeInteractionSameSublattice::calcConstantValues(const Cell &cell)
{
  // cout << "cell check(calcConstantValues): " << cell.begin()->getName() << endl;
  r = cell.getPosition(sl_r);
  M = cell.size();

  // cout << "r: " << r << " " << ", M: " << M << endl;
  double Sr = cell.getSublattice(sl_r).getMoment();

  // cout << "Sr: " << Sr << endl;

  neighbors.findNeighbors(cell, sl_r, sl_r, min, max);
  auto z_rs = static_cast<double>(neighbors.size());
  // cout << "z_rs = " << z_rs << endl;

  // cout << "cell check(calcConstantValues): " << cell2.begin()->getName() << endl;

  LNrr = 0.5 * z_rs * value * Sr;
  LNrs = -0.5 * z_rs * value * Sr;
  // cout << LNrr << " " << LNrs << endl;
}

void ExchangeInteractionSameSublattice::calculateEnergy(const Cell &cell, double &energy)
{
  r = cell.getPosition(sl_r);
  neighbors.findNeighbors(cell, sl_r, sl_r, min, max);
  auto z_rs = static_cast<double>(neighbors.size());
  double Sr = cell[r].getMoment();
  size_t numberOfAtoms = cell[r].size();
  energy -= value * z_rs * numberOfAtoms * Sr * Sr;
}

void ExchangeInteractionSameSublattice::calculateFirstOrderTerms(const Cell & /*cell*/, VectorXcd & /*elements*/)
{
  // first order terms are 0.0
}

void ExchangeInteractionSameSublattice::updateMatrix(const Eigen::Vector3d &K, MatrixXcd &LN) const
{
  complex<double> gamma_rs = neighbors.getGamma(K);
  // cout << value << " " << sl_r << " " << sl_r << " " << gamma_rs << endl;

  // cout << "number of neighbors: " << neighbors.size() << endl;
  // cout << gamma_rs <<  " " << conj(gamma_rs) << endl;
  // cout << LNrr/value << endl;

  LN(r, r) += LNrr + LNrs * conj(gamma_rs);
  LN(r + M, r + M) += LNrr + conj(LNrs) * conj(gamma_rs);
}
}
