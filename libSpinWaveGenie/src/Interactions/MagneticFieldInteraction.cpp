//
//  MagneticFieldInteraction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/17/14.
//
//

#include "SpinWaveGenie/Interactions/MagneticFieldInteraction.h"

using namespace std;

namespace SpinWaveGenie
{

MagneticFieldInteraction::MagneticFieldInteraction(string name_in, double value_in, Vector3 unitVectorIn,
                                                   string sl_r_in)
    : name(name_in), r(0), M(0)
{
  this->updateInteraction(value_in, unitVectorIn, sl_r_in);
}

Interaction *MagneticFieldInteraction::do_clone() const { return new MagneticFieldInteraction(*this); }

void MagneticFieldInteraction::updateInteraction(double value_in, Vector3 unitVectorIn, string sl_r_in)
{
  value = value_in;
  directions = unitVectorIn;
  directions.normalize();
  sl_r = sl_r_in;
}

const string &MagneticFieldInteraction::getName() { return name; }

void MagneticFieldInteraction::updateValue(double value_in) { value = value_in; }

vector<string> MagneticFieldInteraction::sublattices() const
{
  vector<string> sl;
  sl.push_back(sl_r);
  return sl;
}

void MagneticFieldInteraction::calculateEnergy(Cell &cell, double &energy)
{
  if (r < 0)
  {
    r = cell.getPosition(sl_r);
  }
  double S = cell[r].getMoment();
  const Matrix3 &inv = cell[r].getInverseMatrix();

  for (int i = 0; i < 3; i++)
  {
    if (std::abs(directions(i)) > 1.0e-10)
    {
      energy -= value * S * directions(i) * inv(i, 2);
    }
  }
}

void MagneticFieldInteraction::calculateFirstOrderTerms(Cell &cell, Eigen::VectorXcd &elements)
{
  complex<double> XI(0.0, 1.0);
  r = cell.getPosition(sl_r);
  M = cell.size();
  double S = cell[r].getMoment();
  const Matrix3 &inv = cell[r].getInverseMatrix();

  for (int i = 0; i < 3; i++)
  {
    // cout << i << " " << j << " " << directions(i,j) << endl;
    if (abs(directions(i)) > 1.0e-10)
    {
      double X = value * directions(i) * sqrt(S / 2.0);
      complex<double> nu = inv(i, 0) + XI * inv(i, 1);
      // cout << "nu= " << nu << endl;
      elements(r) -= X * conj(nu);
      elements(r + M) -= X * nu;
    }
  }
}

void MagneticFieldInteraction::calcConstantValues(Cell &cell)
{
  complex<double> XI(0.0, 1.0);
  // find location of r
  r = cell.getPosition(sl_r);
  M = cell.size();

  const Matrix3 &inv = cell[r].getInverseMatrix();

  LNrr = complex<double>(0.0, 0.0);

  for (int i = 0; i < 3; i++)
  {
    // cout << i << " " << j << " " << directions(i,j) << endl;
    if (abs(directions(i)) > 1.0e-10)
    {
      LNrr += 0.5 * value * directions(i) * inv(i, 2);
    }
  }
  // cout << LNrr << endl;
}

void MagneticFieldInteraction::updateMatrix(Vector3 K, Eigen::MatrixXcd &LN)
{
  //unused parameter K;
  (void)K;
  LN(r, r) += LNrr;
  LN(r + M, r + M) += LNrr;
}
}
