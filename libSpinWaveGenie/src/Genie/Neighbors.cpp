#include <iostream>
#include "boost/format.hpp"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Containers/Cell.h"

using std::string;
using std::complex;
using std::cout;
using std::endl;

namespace SpinWaveGenie
{
bool Neighbors::empty() { return neighborList.empty(); }

void Neighbors::findNeighbors(const Cell &cell, const string &sl1, const string &sl2, double min, double max)
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
        Vector3 atomdistance(atom2.get<0>() - atom1.get<0>(), atom2.get<1>() - atom1.get<1>(),
                             atom2.get<2>() - atom1.get<2>());
        // cout << "atom2:  " << atom2->get<0>() << " " << atom2->get<1>() << " " << atom2->get<2>() << endl;
        // cout << "atom1:  " << atom1->get<0>() << " " << atom1->get<1>() << " " << atom1->get<2>() << endl;
        // cout << atomdistance.transpose() << endl;
        for (long n1 = -supercellSize; n1 <= supercellSize; n1++)
        {
          for (long n2 = -supercellSize; n2 <= supercellSize; n2++)
          {
            for (long n3 = -supercellSize; n3 <= supercellSize; n3++)
            {
              // find distance between supercells
              Vector3 dispRLU(n1, n2, n3);
              Vector3 dispAng = dispRLU.transpose() * cell.getBasisVectors();
              Vector3 relativeDistance = atomdistance + dispAng;
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
        cout << "Old number of neighbors:" << numberNeighbors << "\n";
        cout << "Old positions:\n";

        for (const auto &MyPosition : neighborList)
        {
          cout << "\t" << MyPosition.get<0>() << " " << MyPosition.get<1>() << " " << MyPosition.get<2>() << "\n";
        }

        cout << "New number of neighbors:" << Neighbors.size() << "\n";
        for (const auto &MyPosition : Neighbors)
        {
          cout << "\t" << MyPosition.get<0>() << " " << MyPosition.get<1>() << " " << MyPosition.get<2>() << "\n";
        }
      }
    }
  }
}

complex<double> Neighbors::getGamma(const Vector3 &K) const
{
  complex<double> MXI(0.0, -1.0);
  complex<double> gamma_rs = complex<double>(0.0, 0.0);
  for (const auto &nbr : neighborList)
  {
    // cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << endl;
    double dot_prod = K[0] * nbr.get<0>() + K[1] * nbr.get<1>() + K[2] * nbr.get<2>();
    gamma_rs += exp(MXI * dot_prod);
    // cout << "gamma_rs = " << gamma_rs << " " << numberNeighbors << endl;
  }
  return gamma_rs / static_cast<double>(numberNeighbors);
}

std::size_t Neighbors::size() { return numberNeighbors; }

Neighbors::Iterator Neighbors::begin() { return neighborList.begin(); }

Neighbors::Iterator Neighbors::end() { return neighborList.end(); }

Neighbors::ConstIterator Neighbors::begin() const { return neighborList.cbegin(); }

Neighbors::ConstIterator Neighbors::end() const { return neighborList.cend(); }

Neighbors::ConstIterator Neighbors::cbegin() const { return neighborList.cbegin(); }

Neighbors::ConstIterator Neighbors::cend() const { return neighborList.cend(); }

std::ostream &operator<<(std::ostream &output, const Neighbors &n)
{
  output << "  x         y         z         r\n";
  for (const auto &nbr : n)
  {
    double dist = sqrt(pow(nbr.get<0>(), 2) + pow(nbr.get<1>(), 2) + pow(nbr.get<2>(), 2));
    output << boost::format("%9.5f %9.5f %9.5f %9.5f\n") % nbr.get<0>() % nbr.get<1>() % nbr.get<2>() % dist;
  }
  return output;
}
}
