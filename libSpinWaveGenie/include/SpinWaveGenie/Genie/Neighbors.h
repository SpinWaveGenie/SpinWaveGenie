#ifndef __Neighbors_H__
#define __Neighbors_H__

#include "Eigen/Core"
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"
#include "SpinWaveGenie/Export.h"
#include <ostream>

namespace SpinWaveGenie
{
class Sublattice;
class Cell;

//! Finds neighbors of two sublattices between distances min and max.
/*!
 The Neighbors class find the relative distance separating neighboring atoms.
 The specific neighbors are described by the two
 sublattices they belong to and the minimum and maximum distance they are separated by.
 */

class SPINWAVEGENIE_EXPORT Neighbors
{
public:
  //! Returns whether of not neighbors have been calculated previously;
  //! \param
  bool empty();
  //! Finds neighbors of two sublattices between distances min and max.
  //! \param cell unit cell
  //! \param sl1 Name of first sublattice
  //! \param sl2 Name of second sublattice
  //! \param min Minimum distance considered, in Angstroms
  //! \param max Maximum distance considered, in Angstroms
  void findNeighbors(const Cell &cell, const std::string &sl1, const std::string &sl2, double min, double max);
  //! Get the number of neighbors.
  std::size_t size() const { return numberNeighbors; }
  //! Get variable \f$ \Gamma = \frac{1}{z_{rs}} \sum_{d} e^{-i \boldmath{k} \cdot \boldmath{d}} \f$
  //! described in J Phys. Condens. Matter 21 216001 (2009)
  //! \param K k vector used in spin wave calculation.
  std::complex<double> getGamma(const Eigen::Vector3d &K) const;
  using Iterator = UniqueThreeVectors<double>::Iterator;
  using ConstIterator = UniqueThreeVectors<double>::ConstIterator;
  //! \return Returns an Iterator pointing to the first element of the neighbor list
  Iterator begin() { return neighborList.begin(); }
  //! \return Returns an Iterator pointing to the final element of the neighbor list
  Iterator end() { return neighborList.end(); }
  //! \return Returns an Iterator pointing to the first element of the neighbor list
  ConstIterator begin() const { return neighborList.cbegin(); }
  //! \return Returns an Iterator pointing to the final element of the neighbor list
  ConstIterator end() const { return neighborList.cend(); }
  //! \return Returns an ConstIterator pointing to the first element of the neighbor list
  ConstIterator cbegin() const { return neighborList.cbegin(); }
  //! \return Returns an ConstIterator pointing to the final element of the neighbor list
  ConstIterator cend() const { return neighborList.cend(); }
  friend SPINWAVEGENIE_EXPORT std::ostream &operator<<(std::ostream &output, const Neighbors &n);

private:
  UniqueThreeVectors<double> neighborList;
  std::size_t numberNeighbors{0};
};
} // namespace SpinWaveGenie
#endif // __Neighbors_H__
