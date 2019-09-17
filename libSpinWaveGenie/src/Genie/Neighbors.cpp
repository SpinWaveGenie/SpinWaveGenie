#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Containers/Cell.h"
#include "boost/format.hpp"
#include <iostream>

namespace SpinWaveGenie
{
bool Neighbors::empty() { return neighborList.empty(); }

void Neighbors::findNeighbors(const Cell &cell, const std::string &sl1, const std::string &sl2, double min, double max)
{
  // In principle, we only need to iterate over one atom in the first sublattice. However, iterating over
  // all atoms provides a good check that all atoms have the same number of neighbors in the same relative
  // positions

  bool firstTime = true;
  for (const auto &atom1 : cell.getSublattice(sl1))
  // const auto &atom1 = *(cell.getSublattice(sl1).begin());
  {
    //  A 5x5x5 supercell should good enough for any physical interaction.
    UniqueThreeVectors<double> Neighbors;
    for (long supercellSize = 5; supercellSize <= 5; supercellSize++)
    {
      // cout << supercellSize << endl;
      for (const auto &atom2 : cell.getSublattice(sl2))
      {
        Eigen::Vector3d atomdistance(atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]);
        // cout << "atom2:  " << atom2->get<0>() << " " << atom2->get<1>() << " " << atom2->get<2>() << endl;
        // cout << "atom1:  " << atom1->get<0>() << " " << atom1->get<1>() << " " << atom1->get<2>() << endl;
        // cout << atomdistance.transpose() << endl;
        for (int32_t n1 = -supercellSize; n1 <= supercellSize; n1++)
        {
          for (int32_t n2 = -supercellSize; n2 <= supercellSize; n2++)
          {
            for (int32_t n3 = -supercellSize; n3 <= supercellSize; n3++)
            {
              // find distance between supercells
              Eigen::Vector3d dispRLU(n1, n2, n3);
              Eigen::Vector3d dispAng = dispRLU.transpose() * cell.getBasisVectors();
              Eigen::Vector3d relativeDistance = atomdistance + dispAng;
              double norm = relativeDistance.norm();
              // check if norm is between min and max
              if (norm < max && norm > min)
              {
                // cout << "candidate: " << relativeDistance.transpose() << endl;
                Neighbors.insert(relativeDistance[0], relativeDistance[1], relativeDistance[2]);
              }
            }
          }
        }
      }
    }
    if (firstTime)
    {
      firstTime = false;
      numberNeighbors = Neighbors.size();
      neighborList = Neighbors;
    }
    else
    {
      if (!(Neighbors == neighborList))
      {
        std::cout << "Old number of neighbors:" << numberNeighbors << "\n";
        std::cout << "Old positions:\n";

        for (const auto &MyPosition : neighborList)
        {
          std::cout << "\t" << MyPosition[0] << " " << MyPosition[1] << " " << MyPosition[2] << "\n";
        }

        std::cout << "New number of neighbors:" << Neighbors.size() << "\n";
        for (const auto &MyPosition : Neighbors)
        {
          std::cout << "\t" << MyPosition[0] << " " << MyPosition[1] << " " << MyPosition[2] << "\n";
        }
      }
    }
  }
}

std::complex<double> Neighbors::getGamma(const Eigen::Vector3d &K) const
{
  std::complex<double> MXI(0.0, -1.0);
  std::complex<double> gamma_rs = std::complex<double>(0.0, 0.0);
  for (const auto &nbr : neighborList)
  {
    // cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << endl;
    double dot_prod = K[0] * nbr[0] + K[1] * nbr[1] + K[2] * nbr[2];
    gamma_rs += exp(MXI * dot_prod);
    // cout << "gamma_rs = " << gamma_rs << " " << numberNeighbors << endl;
  }
  return gamma_rs / static_cast<double>(numberNeighbors);
}

std::ostream &operator<<(std::ostream &output, const Neighbors &n)
{
  output << "  x         y         z         r\n";
  for (const auto &nbr : n)
  {
    double dist = sqrt(pow(nbr[0], 2) + pow(nbr[1], 2) + pow(nbr[2], 2));
    output << boost::format("%9.5f %9.5f %9.5f %9.5f\n") % nbr[0] % nbr[1] % nbr[2] % dist;
  }
  return output;
}
} // namespace SpinWaveGenie
