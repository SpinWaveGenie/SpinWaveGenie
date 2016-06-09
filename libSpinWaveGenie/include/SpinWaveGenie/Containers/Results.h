//
//  Results.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/1/14.
//
//

#ifndef __spin_wave_genie__Results__
#define __spin_wave_genie__Results__

#include "SpinWaveGenie/Export.h"
#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

//! Stores a frequency and intensity pair.
/*!
 This struct contains a calculated frequency and its associated intensity.
 The operator< is used to sort the Points by frequency.
 */
struct Point
{
  //! Frequency associated with a given excitation (in meV).
  double frequency;
  //! Measurable intensity of a given excitation (arb. units).
  double intensity;
  //! \return Whether this Point has lower frequency than Point val.
  bool operator<(const Point &val) const { return frequency < val.frequency; }
};

//! Stores results of SpinWaveGenie Calculations at the given Q-point.
/*!
 This container stores the results of a given SpinWaveGenie calculation. Each "Point"
 contains a calculated frequency and intensity. Results can be sorted by frequency,
 or filtered so that branches with identical frequencies are combined, or branches without
 significant intensity are discarded.
 */

class SPINWAVEGENIE_EXPORT Results
{
public:
  //! Insert Point struct into container.
  void insert(Point value);
  //! \return size of Results container.
  std::size_t size() const { return results.size(); }
  //! Sort Results by frequency.
  void sort();
  //! Clear results container so that the size is zero.
  void clear();
  //! Filter Points by combining those with identical frequencies.
  void uniqueSolutions();
  //! Filter Points by removing those without significant intensity.
  void significantSolutions(double ETS = 1.0e-10);
  const Point &operator[](std::size_t bin) { return results[bin]; }
  typedef std::vector<Point>::iterator Iterator;
  typedef std::vector<Point>::const_iterator ConstIterator;
  //! \return Returns an Iterator pointing to the first element in the container.
  Iterator begin() { return results.begin(); }
  //! \return Returns an Iterator pointing to the end of the container.
  Iterator end() { return results.end(); }
  //! \return Returns an ConstIterator pointing to the first element in the container.
  ConstIterator begin() const { return results.cbegin(); }
  //! \return Returns an ConstIterator pointing to the end of the container.
  ConstIterator end() const { return results.cend(); }
  //! \return Returns an ConstIterator pointing to the first element in the container.
  ConstIterator cbegin() const { return results.cbegin(); }
  //! \return Returns an ConstIterator pointing to the end of the container.
  ConstIterator cend() const { return results.cend(); }
  friend SPINWAVEGENIE_EXPORT std::ostream &operator<<(std::ostream &output, const Results &n);

private:
  std::vector<Point> results;
};
}
#endif /* defined(__spin_wave_genie__Results__) */
