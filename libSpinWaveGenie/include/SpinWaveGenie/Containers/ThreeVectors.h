#ifndef __ThreeVectors__
#define __ThreeVectors__

#include <array>
#include <cassert>
#include <iostream>
#include <vector>

//! Array of Struct used for storing vectors with three components.
/*!
 Vectors containing three elements */

namespace SpinWaveGenie
{

template <typename T> class ThreeVectors
{

public:
  bool empty() const;
  //! insert three elements x,y,z
  //! \param x zeroth element of type T
  //! \param y first element of type T
  //! \param z second element of type T
  bool insert(T x, T y, T z);
  using size_type = typename std::vector<std::array<T, 3>>::size_type;
  using difference_type = typename std::vector<std::array<T, 3>>::difference_type;
  using Iterator = typename std::vector<std::array<T, 3>>::iterator;
  using ConstIterator = typename std::vector<std::array<T, 3>>::const_iterator;
  //! \return number of elements in the ThreeVector
  size_type size() const;
  //! Clears all data stored in the ThreeVector.
  void clear();
  //! \return Returns an Iterator pointing to the first element
  Iterator begin();
  //! \return Returns an Iterator pointing to the end of the vector
  Iterator end();
  //! \return Returns a ConstIterator pointing to the first element
  ConstIterator begin() const;
  //! \return Returns an ConstIterator pointing to the end of the vector
  ConstIterator end() const;
  //! \return Returns a ConstIterator pointing to the first element
  ConstIterator cbegin() const;
  //! \return Returns an ConstIterator pointing to the end of the vector
  ConstIterator cend() const;

protected:
  std::vector<std::array<T, 3>> values;
};

template <typename T> bool ThreeVectors<T>::empty() const { return values.empty(); }

template <typename T> bool ThreeVectors<T>::insert(T x, T y, T z)
{
  values.push_back({{x,y,z}});
  return true;
}

template <typename T> typename ThreeVectors<T>::size_type ThreeVectors<T>::size() const { return values.size(); }

template <typename T> typename ThreeVectors<T>::Iterator ThreeVectors<T>::begin()
{
  return values.begin();
}

template <typename T> typename ThreeVectors<T>::Iterator ThreeVectors<T>::end()
{
  return values.end();
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::begin() const
{
  return values.cbegin();
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::end() const
{
  return values.cend();
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cbegin() const
{
  return values.cbegin();
}

template <typename T> typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cend() const
{
  return values.cend();
}

template <typename T> void ThreeVectors<T>::clear()
{
  values.clear();
}
}
#endif /* defined(__ThreeVectors__) */
