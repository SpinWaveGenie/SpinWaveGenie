#ifndef __ThreeVectors__
#define __ThreeVectors__

#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

template<class T>
class ThreeVectors
{
protected:
    typedef typename std::vector<T>::iterator ValueIterator;
    typedef typename std::vector<T>::const_iterator ConstValueIterator;
public:
    bool operator==(ThreeVectors& other);
    void insert(T x, T y,T z);
    typedef boost::zip_iterator<boost::tuple<ValueIterator,ValueIterator,ValueIterator> > Iterator;
    typedef boost::zip_iterator<boost::tuple<ConstValueIterator,ConstValueIterator,ConstValueIterator> > ConstIterator;
    size_t size();
    void clear();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator begin();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator end();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    ConstIterator cbegin();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    ConstIterator cend();
protected:
    std::vector<T> valuesX;
    std::vector<T> valuesY;
    std::vector<T> valuesZ;
};

template<class T>
bool ThreeVectors<T>::operator==(ThreeVectors<T>& other)
{
    for(auto MyThreeVectors = this->begin(); MyThreeVectors!= this->end(); MyThreeVectors++)
    {
        bool ContainsPosition = false;
        for(auto OtherThreeVectors = other.begin(); OtherThreeVectors!= other.end(); OtherThreeVectors++)
        {
            double epsilon = 1.0e-8;
            T x = MyThreeVectors->template get<0>() - OtherThreeVectors->template get<0>();
            T y = MyThreeVectors->template get<1>() - OtherThreeVectors->template get<1>();
            T z = MyThreeVectors->template get<2>() - OtherThreeVectors->template get<2>();
            if (std::abs(x) < epsilon && std::abs(y) < epsilon && std::abs(z) < epsilon)
                ContainsPosition = true;
        }
        if(ContainsPosition == false)
            return ContainsPosition;
    }
    return true;
}

template<class T>
void ThreeVectors<T>::insert(T x,T y,T z)
{
    valuesX.push_back(x);
    valuesY.push_back(y);
    valuesZ.push_back(z);
}

template<class T>

size_t ThreeVectors<T>::size()
{
    return valuesX.size();
}

template<class T>
typename ThreeVectors<T>::Iterator ThreeVectors<T>::begin()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.begin(),valuesY.begin(),valuesZ.begin()));
}

template<class T>
typename ThreeVectors<T>::Iterator ThreeVectors<T>::end()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.end(),valuesY.end(),valuesZ.end()));
}

template<class T>
typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cbegin()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.cbegin(),valuesY.cbegin(),valuesZ.cbegin()));
}

template<class T>
typename ThreeVectors<T>::ConstIterator ThreeVectors<T>::cend()
{
    return boost::make_zip_iterator(boost::make_tuple(valuesX.cend(),valuesY.cend(),valuesZ.cend()));
}

template<class T>
void ThreeVectors<T>::clear()
{
    valuesX.clear();
    valuesY.clear();
    valuesZ.clear();
}

#endif /* defined(__ThreeVectors__) */
