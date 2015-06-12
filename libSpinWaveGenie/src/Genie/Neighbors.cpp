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
bool Neighbors::empty() { return neighborList.empty(); };

void Neighbors::findNeighbors(Cell &cell, string sl1, string sl2, double min, double max)
{
  // In principle, we only need to iterate over one atom in the first sublattice. However, iterating over
  // all atoms provides a good check that all atoms have the same number of neighbors in the same relative
  // positions

  for (Sublattice::Iterator atom1 = cell.getSublattice(sl1).begin(); atom1 != cell.getSublattice(sl1).end(); ++atom1)
  // auto atom1 = cell.getSublattice(sl1).begin();
  {
    //  A 5x5x5 supercell should good enough for any physical interaction.
    UniqueThreeVectors<double> Neighbors;
    for (long supercellSize = 5; supercellSize <= 5; supercellSize++)
    {
      // cout << supercellSize << endl;
      for (Sublattice::Iterator atom2 = cell.getSublattice(sl2).begin(); atom2 != cell.getSublattice(sl2).end();
           ++atom2)
      {
        Vector3 atomdistance(atom2->get<0>() - atom1->get<0>(), atom2->get<1>() - atom1->get<1>(),
                             atom2->get<2>() - atom1->get<2>());
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
    if (atom1 == cell.getSublattice(sl1).begin())
    {
      numberNeighbors = Neighbors.size();
      neighborList = Neighbors;
    }
    else
    {
      if (!(Neighbors == neighborList))
      {
        cout << "Old number of neighbors:" << numberNeighbors << endl;
        cout << "Old positions:" << endl;
        for (auto MyPosition = neighborList.begin(); MyPosition != neighborList.end(); MyPosition++)
        {
          cout << "\t" << MyPosition->get<0>() << " " << MyPosition->get<1>() << " " << MyPosition->get<2>() << endl;
        }
        cout << "New number of neighbors:" << Neighbors.size() << endl;
        for (auto MyPosition = Neighbors.begin(); MyPosition != Neighbors.end(); MyPosition++)
        {
          cout << "\t" << MyPosition->get<0>() << " " << MyPosition->get<1>() << " " << MyPosition->get<2>() << endl;
        }
      }
    }
  }
}

complex<double> Neighbors::getGamma(Vector3 K)
{
  complex<double> MXI(0.0, -1.0);
  complex<double> gamma_rs = complex<double>(0.0, 0.0);
  for (Iterator nbr = this->begin(); nbr != this->end(); ++nbr)
  {
    // cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << endl;
    double dot_prod = K[0] * nbr->get<0>() + K[1] * nbr->get<1>() + K[2] * nbr->get<2>();
    gamma_rs += exp(MXI * dot_prod);
    // cout << "gamma_rs = " << gamma_rs << " " << numberNeighbors << endl;
  }
  return gamma_rs / static_cast<double>(numberNeighbors);
}

std::size_t Neighbors::size() { return numberNeighbors; }

Neighbors::Iterator Neighbors::begin() { return Iterator(neighborList.begin()); }

Neighbors::Iterator Neighbors::end() { return Iterator(neighborList.end()); }

Neighbors::ConstIterator Neighbors::cbegin() const { return ConstIterator(neighborList.cbegin()); }

Neighbors::ConstIterator Neighbors::cend() const { return ConstIterator(neighborList.cend()); }

std::ostream &operator<<(std::ostream &output, const Neighbors &n)
{
  output << "  x         y         z\n";
  for (auto nbr = n.cbegin(); nbr != n.cend(); ++nbr)
  {
    output << boost::format("%9.5f %9.5f %9.5f\n") % nbr->get<0>() % nbr->get<1>() % nbr->get<2>();
  }
  return output;
}
}
