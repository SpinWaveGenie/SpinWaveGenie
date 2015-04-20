#ifndef __Neighbors_H__
#define __Neighbors_H__

#define _USE_MATH_DEFINES
#include <ostream>
#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"

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

class Neighbors
{
public:
  Neighbors() : numberNeighbors(0){};
  //! Returns whether of not neighbors have been calculated previously;
  //! \param
  bool empty();
  //! Finds neighbors of two sublattices between distances min and max.
  //! \param cell unit cell
  //! \param sl1 Name of first sublattice
  //! \param sl2 Name of second sublattice
  //! \param min Minimum distance considered, in Angstroms
  //! \param max Maximum distance considered, in Angstroms
  void findNeighbors(Cell &cell, std::string sl1, std::string sl2, double min, double max);
  //! Get the number of neighbors.
  double size();
  //! Get variable \f$ \Gamma = \frac{1}{z_{rs}} \sum_{d} e^{-i \boldmath{k} \cdot \boldmath{d}} \f$
  //! described in J Phys. Condens. Matter 21 216001 (2009)
  //! \param K k vector used in spin wave calculation.
  std::complex<double> getGamma(Vector3 K);
  typedef UniqueThreeVectors<double>::Iterator Iterator;
  typedef UniqueThreeVectors<double>::ConstIterator ConstIterator;
  //! \return Returns an Iterator pointing to the first element of the neighbor list
  Iterator begin();
  //! \return Returns an Iterator pointing to the final element of the neighbor list
  Iterator end();
  //! \return Returns an ConstIterator pointing to the first element of the neighbor list
  ConstIterator cbegin() const;
  //! \return Returns an ConstIterator pointing to the final element of the neighbor list
  ConstIterator cend() const;
  friend std::ostream& operator<<( std::ostream &output, const Neighbors &n);
private:
  UniqueThreeVectors<double> neighborList;
  double numberNeighbors;
};
}
#endif // __Neighbors_H__
