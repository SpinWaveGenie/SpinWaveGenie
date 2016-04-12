#include <iostream>
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Interactions/ExchangeInteraction.h"
#include "SpinWaveGenie/Memory.h"

using namespace std;
using namespace Eigen;

namespace SpinWaveGenie
{

ExchangeInteraction::ExchangeInteraction(string name_in, double value_in, string sl_r_in, string sl_s_in, double min_in,
                                         double max_in)
    : name(name_in), r(0), s(0), M(0)
{
  name = name_in;
  this->updateInteraction(value_in, sl_r_in, sl_s_in, min_in, max_in);
}

std::unique_ptr<Interaction> ExchangeInteraction::clone() const
{
  return memory::make_unique<ExchangeInteraction>(*this);
}

void ExchangeInteraction::updateInteraction(double value_in, string sl_r_in, string sl_s_in, double min_in,
                                            double max_in)
{
  value = value_in;
  sl_r = sl_r_in;
  sl_s = sl_s_in;
  min = min_in;
  max = max_in;
}

const string &ExchangeInteraction::getName() { return name; }

void ExchangeInteraction::updateValue(double value_in) { value = value_in; }

std::array<std::string, 2> ExchangeInteraction::sublattices() const
{
  if (sl_r < sl_s)
    return {{sl_r, sl_s}};
  else
    return {{sl_s, sl_r}};
}

void ExchangeInteraction::calcConstantValues(Cell &cell)
{
  r = cell.getPosition(sl_r);
  s = cell.getPosition(sl_s);
  M = cell.size();

  double Sr = cell.getSublattice(sl_r).getMoment();
  double Ss = cell.getSublattice(sl_s).getMoment();

  Matrix3 Frs = cell.getSublattice(sl_r).getRotationMatrix() * cell.getSublattice(sl_s).getInverseMatrix();

  Matrix3 Fsr = cell.getSublattice(sl_s).getRotationMatrix() * cell.getSublattice(sl_r).getInverseMatrix();

  neighbors.findNeighbors(cell, sl_r, sl_s, min, max);
  double z_rs = static_cast<double>(neighbors.size());

  complex<double> G1rs = -0.5 * complex<double>(Frs(0, 0) + Frs(1, 1), Frs(1, 0) - Frs(0, 1));
  complex<double> G2rs = -0.5 * complex<double>(Frs(0, 0) - Frs(1, 1), -Frs(1, 0) - Frs(0, 1));
  complex<double> G1sr = -0.5 * complex<double>(Fsr(0, 0) + Fsr(1, 1), Fsr(1, 0) - Fsr(0, 1));
  complex<double> G2sr = -0.5 * complex<double>(Fsr(0, 0) - Fsr(1, 1), -Fsr(1, 0) - Fsr(0, 1));

  LNrr = 0.25 * z_rs * value * Ss * (Frs(2, 2) + Fsr(2, 2));
  LNss = 0.25 * z_rs * value * Sr * (Frs(2, 2) + Fsr(2, 2));
  LNrs = 0.25 * z_rs * value * sqrt(Sr * Ss) * (G1rs + conj(G1sr));
  LNrsM = 0.25 * z_rs * value * sqrt(Sr * Ss) * (conj(G2rs) + conj(G2sr));

  // cout << LNrr << " "<< LNss << " " << LNrs << " " << LNrsM << endl;
}

void ExchangeInteraction::calculateEnergy(Cell &cell, double &energy)
{
  if (neighbors.empty())
  {
    neighbors.findNeighbors(cell, sl_r, sl_s, min, max);
    r = cell.getPosition(sl_r);
    s = cell.getPosition(sl_s);
  }
  double z_rs = static_cast<double>(neighbors.size());

  double Sr = cell[r].getMoment();
  double Ss = cell[s].getMoment();

  Matrix3 Frs = cell[r].getRotationMatrix() * cell[s].getInverseMatrix();

  Matrix3 Fsr = cell[s].getRotationMatrix() * cell[r].getInverseMatrix();

  energy -= 0.5 * value * z_rs * Sr * Ss * (Frs(2, 2) + Fsr(2, 2));
}

void ExchangeInteraction::calculateFirstOrderTerms(Cell &cell, VectorXcd &elements)
{
  r = cell.getPosition(sl_r);
  s = cell.getPosition(sl_s);
  M = cell.size();

  double Sr = cell[r].getMoment();
  double Ss = cell[s].getMoment();

  neighbors.findNeighbors(cell, sl_r, sl_s, min, max);
  double z_rs = static_cast<double>(neighbors.size());

  Matrix3 Frs = cell[r].getRotationMatrix() * cell[s].getInverseMatrix();

  Matrix3 Fsr = cell[s].getRotationMatrix() * cell[r].getInverseMatrix();

  complex<double> F1rs(Frs(0, 2), Frs(1, 2));
  complex<double> F2rs(Frs(2, 0), Frs(2, 1));
  complex<double> F1sr(Fsr(0, 2), Fsr(1, 2));
  complex<double> F2sr(Fsr(2, 0), Fsr(2, 1));

  complex<double> LNr = -1.0 * sqrt(Sr) * Ss / (2.0 * sqrt(2.0)) * z_rs * value * conj(F1rs + F2sr);
  complex<double> LNs = -1.0 * sqrt(Ss) * Sr / (2.0 * sqrt(2.0)) * z_rs * value * conj(F1sr + F2rs);

  elements[r] += LNr;
  elements[s] += LNs;
  elements[r + M] += conj(LNr);
  elements[s + M] += conj(LNs);
}

void ExchangeInteraction::updateMatrix(Vector3d K, MatrixXcd &LN) const
{
  std::complex<double> gamma_rs = neighbors.getGamma(K);
  // cout << sl_r << " " << sl_s << " " << gamma_rs << endl;

  LN(r, r) += LNrr;
  LN(r, s) += LNrs * conj(gamma_rs);
  LN(r, s + M) += LNrsM * conj(gamma_rs);

  LN(s, r) += conj(LNrs) * gamma_rs;
  LN(s, s) += LNss;
  LN(s, r + M) += LNrsM * gamma_rs;

  LN(r + M, s) += conj(LNrsM) * conj(gamma_rs);
  LN(r + M, r + M) += LNrr;
  LN(r + M, s + M) += conj(LNrs) * conj(gamma_rs);

  LN(s + M, r) += conj(LNrsM) * gamma_rs;
  LN(s + M, r + M) += LNrs * gamma_rs;
  LN(s + M, s + M) += LNss;
}
}
