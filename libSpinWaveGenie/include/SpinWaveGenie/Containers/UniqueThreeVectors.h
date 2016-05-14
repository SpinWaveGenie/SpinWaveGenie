//
//  UniquePositions.h
//  positions
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#ifndef __positions__UniqueThreeVectors__
#define __positions__UniqueThreeVectors__

#include <complex>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/range/iterator_range.hpp>
#include "SpinWaveGenie/Containers/ThreeVectors.h"

namespace SpinWaveGenie
{

//! Structure of Arrays used for storing vectors with three components. Only unique vectors are stored.
/*!
 Vectors containing three elements are not ideally coellesced in memory.
 Therefore, we use an alternative where where each component is stored
 in a separate vector and accessed using an Iterator. Unlike ThreeVectors, only the
 unique vectors are stored.
 */

template <typename T> class UniqueThreeVectors : public ThreeVectors<T>
{
public:
  //!
  bool operator==(const UniqueThreeVectors &other);
  bool insert(T x, T y, T z);
};

template <typename T> bool UniqueThreeVectors<T>::insert(T x, T y, T z)
{
  std::array<T, 3> meow = {{x, y, z}};
  if (std::find(this->begin(), this->end(), meow) == this->end())
  {
    ThreeVectors<T>::insert(x, y, z);
    return true;
  }
  else
  {
    return false;
  }
}

template <typename T> bool UniqueThreeVectors<T>::operator==(const UniqueThreeVectors<T> &other)
{
  for (const auto &value : *this)
  {
    if (std::find(other.cbegin(), other.cend(), value) == other.cend())
      return false;
  }
  return true;
}
}
#endif /* defined(__positions__UniquePositions__) */
