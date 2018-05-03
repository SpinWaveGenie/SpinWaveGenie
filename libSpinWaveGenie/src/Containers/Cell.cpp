#include "SpinWaveGenie/Containers/Cell.h"

#include "boost/math/special_functions/pow.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdexcept>

namespace SpinWaveGenie
{

struct CompareSublatticeNames
{
  CompareSublatticeNames(const std::string &sublatticeName) : name(sublatticeName) {}
  bool operator()(const Sublattice &arg) { return name == arg.getName(); }
  std::string name;
};

void Cell::setBasisVectors(double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg)
{
  //! <a href=https://github.com/mantidproject/documents/blob/master/Design/UBMatriximplementationnotes.pdf> Reference
  //</a>
  constexpr double deg2rad = M_PI / 180.0;
  double alpha, beta, gamma;
  alpha = alpha_deg * deg2rad;
  beta = beta_deg * deg2rad;
  gamma = gamma_deg * deg2rad;

  double ci = c * std::cos(beta);
  double cj = c * (std::cos(alpha) - std::cos(gamma) * std::cos(beta)) / std::sin(gamma);
  double ck = c / std::sin(gamma) *
              std::sqrt(1.0 - boost::math::pow<2>(cos(alpha)) - boost::math::pow<2>(cos(beta)) -
                        boost::math::pow<2>(cos(gamma)) + 2.0 * cos(alpha) * cos(beta) * cos(gamma));

  basisVectors << a, 0.0, 0.0, b *cos(gamma), b * sin(gamma), 0.0, ci, cj, ck;

  // cout << "basis vectors equal" <<basisVectors << endl;

  // basisVectors << a/2.0,a*sqrt(3.0)/2.0,0.0,
  //                a/2.0,-1.0*a*sqrt(3.0)/2.0,0.0,
  //                0.0,0.0,c;

  reciprocalVectors = 2.0 * M_PI * basisVectors.inverse().transpose();

  // cout << "recip vectors equal" <<reciprocalVectors << endl;
}

void Cell::setBasisVectors(double scale, const Eigen::Matrix3d &basis) { basisVectors = scale * basis; }

const Eigen::Matrix3d &Cell::getBasisVectors() const { return basisVectors; }

const Eigen::Matrix3d &Cell::getReciprocalVectors() const { return reciprocalVectors; }

void Cell::addSublattice(const Sublattice &sl)
{
  const std::string &name = sl.getName();
  auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
  if (it != sublatticeInfo.end())
  {
    throw std::invalid_argument("sublattice already defined");
  }
  else
  {
    sublatticeInfo.push_back(sl);
  }
}

Sublattice &Cell::getSublattice(const std::string &name)
{
  auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
  if (it == sublatticeInfo.end())
  {
    throw std::invalid_argument("sublattice not found");
  }
  return *it;
}

const Sublattice &Cell::getSublattice(const std::string &name) const
{
  auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
  if (it == sublatticeInfo.end())
  {
    throw std::invalid_argument("sublattice not found");
  }
  return *it;
}

const Sublattice &Cell::operator[](std::vector<Sublattice>::size_type position) const
{
  return sublatticeInfo[position];
}

std::vector<Sublattice>::difference_type Cell::getPosition(const std::string &name) const
{
  auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
  if (it == sublatticeInfo.end())
  {
    throw std::invalid_argument("sublattice not found");
  }
  return std::distance(sublatticeInfo.begin(), it);
}

void Cell::addAtom(const std::string &name, double x, double y, double z)
{
  Eigen::Vector3d scaled_position(x, y, z);

  // cout << "scaled= " << scaled_position.transpose() << endl;
  // cout << basisVectors << endl;

  Eigen::Vector3d pos = scaled_position.transpose() * basisVectors;

  // cout << "unscaled= " << pos.transpose() << endl;
  // cout << " " <<endl;

  this->getSublattice(name).addAtom(pos[0], pos[1], pos[2]);
}

}
