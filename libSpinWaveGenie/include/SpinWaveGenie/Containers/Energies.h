//
//  Energies.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/27/14.
//
//

#ifndef __spin_wave_genie__Energies__
#define __spin_wave_genie__Energies__

#include "SpinWaveGenie/Export.h"
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

//! Stores values of energy calculated by the SpinWavePlot routines.
/*!
    This container stores values of energy to be calculated by the SpinWavePlot routines.
    Energies are sorted for low to high and for a given energy, one can get the lower
    and upper bound of the index.
 */

class SPINWAVEGENIE_EXPORT Energies
{
public:
  //! Default constructor
  Energies() = default;
  //! Constructs numberPoints evenly spaced energies between minimum and maximum.
  //! \param minimum Minimum energy in range, inclusive.
  //! \param maximum Maximum energy in range, inclusive.
  //! \param numberPoints Number of points in range.
  Energies(double minimum, double maximum, std::size_t numberPoints);
  //! Insert an energy value into container.
  //! \param energy energy value.
  void insert(double energy);
  //! Get index of the first element in the range that is not less than (i.e. greater or equal to) the specified energy.
  //! \param energy energy value.
  //! \return index of lower bound.
  std::size_t getLowerBound(double energy);
  //! Get index of the first element in the range which compares greater than the specified energy.
  //! \param energy energy value.
  //! \return index of upper bound.
  std::size_t getUpperBound(double energy);
  //! Get the number of values stored in this container.
  //! \return size of container.
  std::size_t size() const { return energies.size(); }
  //! Clear all energy values currently stored in this container. this->size() should be zero.
  void clear();
  //! Get C-style pointer to first element in container.
  //! \return pointer to internal vector.
  double *data();
  //! Permit access via the subscript operator.
  //! \param position Index of energy value.
  //! \return Energy value.
  const double &operator[](std::size_t bin) { return energies[bin]; }
  using Iterator = std::vector<double>::iterator;
  using ConstIterator = std::vector<double>::const_iterator;
  //! \return Returns an Iterator pointing to the first element of the container.
  Iterator begin() { return energies.begin(); }
  //! \return Returns an Iterator pointing to the end of the container.
  Iterator end() { return energies.end(); }
  //! \return Returns a ConstIterator pointing to the first element of the container.
  ConstIterator begin() const { return energies.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the container.
  ConstIterator end() const { return energies.cend(); }
  //! \return Returns a ConstIterator pointing to the first element of the container.
  ConstIterator cbegin() const { return energies.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the container.
  ConstIterator cend() const { return energies.cend(); }
  friend SPINWAVEGENIE_EXPORT std::ostream &operator<<(std::ostream &output, const Energies &n);

private:
  std::vector<double> energies;
};

} // namespace SpinWaveGenie
#endif /* defined(__spin_wave_genie__Energies__) */
