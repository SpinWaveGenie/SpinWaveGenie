#ifndef __Positions__
#define __Positions__

#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

class Positions
{
protected:
    typedef std::vector<double>::iterator ValueIterator;
    typedef std::vector<double>::const_iterator ConstValueIterator;
public:
    bool operator==(Positions& other);
    void insert(double x, double y, double z);
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
    std::vector<double> valuesX;
    std::vector<double> valuesY;
    std::vector<double> valuesZ;
};

#endif /* defined(__Positions__) */
