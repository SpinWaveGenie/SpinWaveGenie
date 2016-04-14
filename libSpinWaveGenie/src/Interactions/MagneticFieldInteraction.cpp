//
//  MagneticFieldInteraction.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/17/14.
//
//

#include "SpinWaveGenie/Interactions/MagneticFieldInteraction.h"
#include "SpinWaveGenie/Memory.h"

using namespace std;

namespace SpinWaveGenie
{

MagneticFieldInteraction::MagneticFieldInteraction(string name_in, double value_in, const Vector3 &unitVectorIn,
                                                   const string &sl_r_in)
    : name(std::move(name_in)), r(0), M(0)
{
  this->updateInteraction(value_in, unitVectorIn, sl_r_in);
}

std::unique_ptr<Interaction> MagneticFieldInteraction::clone() const
{
  return memory::make_unique<MagneticFieldInteraction>(*this);
}

void MagneticFieldInteraction::updateInteraction(double value_in, const Vector3 &unitVectorIn, const string &sl_r_in)
{
  value = value_in;
  directions = unitVectorIn;
  directions.normalize();
  sl_r = sl_r_in;
}

const string &MagneticFieldInteraction::getName() const { return name; }

void MagneticFieldInteraction::updateValue(double value_in) { value = value_in; }

std::array<std::string, 2> MagneticFieldInteraction::sublattices() const { return {{sl_r, sl_r}}; }

void MagneticFieldInteraction::calculateEnergy(const Cell &cell, double &energy)
{
  r = cell.getPosition(sl_r);
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

void MagneticFieldInteraction::calculateFirstOrderTerms(const Cell &cell, Eigen::VectorXcd &elements)
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

void MagneticFieldInteraction::calcConstantValues(const Cell &cell)
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

void MagneticFieldInteraction::updateMatrix(const Vector3 & /*K*/, Eigen::MatrixXcd &LN) const
{
  LN(r, r) += LNrr;
  LN(r + M, r + M) += LNrr;
}
}
