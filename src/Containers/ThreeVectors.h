#ifndef __ThreeVectors__
#define __ThreeVectors__

#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

//! Struture of Arrays used for storing vectors with three components.
/*!
 Vectors containing three elements are not ideally coellesced in memory. Therefore, we use an alternative where
 where each component is stored in a separate vector and accessed using an Iterator.
 */

template<typename T>
class ThreeVectors
{
protected:
    typedef typename std::vector<T>::iterator ValueIterator;
    typedef typename std::vector<T>::const_iterator ConstValueIterator;
public:
    //! insert three elements x,y,z
    //! \param x zeroth element of type T
    //! \param y first element of type T
    //! \param z second element of type T
    void insert(T x, T y,T z);
    typedef boost::zip_iterator<boost::tuple<ValueIterator,ValueIterator,ValueIterator> > Iterator;
    typedef boost::zip_iterator<boost::tuple<ConstValueIterator,ConstValueIterator,ConstValueIterator> > ConstIterator;
    //! \return number of elements in the ThreeVector
    size_t size();
    //! Clears all data stored in the ThreeVector.
    void clear();
    //! \return Returns an iterator pointing to the first element
    Iterator begin();
    //! \return Returns an iterator pointing to the end of the vector
    Iterator end();
    //! \return Returns an iterator pointing to the first element
    ConstIterator cbegin();
    //! \return Returns an iterator pointing to the end fo the vector
    ConstIterator cend();
protected:
    std::vector<T> valuesX;
    std::vector<T> valuesY;
    std::vector<T> valuesZ;
};

template<typename T>
void ThreeVectors<T>::insert(T x,T y,T z)
{
    valuesX.push_back(x);
    valuesY.push_back(y);
    valuesZ.push_back(z);
}

template<typename T>

size_t ThreeVectors<T>::size()
{
    return valuesX.size();
}

template<typename T>
typename ThreeVectors<T>::Iterator ThreeVectors<T>::begin()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.begin(),valuesY.begin(),valuesZ.begin()));
}

template<typename T>
typename ThreeVectors<T>::Iterator ThreeVectors<T>::end()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.end(),valuesY.end(),valuesZ.end()));
}

template<typename T>
typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cbegin()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.cbegin(),valuesY.cbegin(),valuesZ.cbegin()));
}

template<typename T>
typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cend()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.cend(),valuesY.cend(),valuesZ.cend()));
}

template<typename T>
void ThreeVectors<T>::clear()
{
    valuesX.clear();
    valuesY.clear();
    valuesZ.clear();
}

#endif /* defined(__ThreeVectors__) */
