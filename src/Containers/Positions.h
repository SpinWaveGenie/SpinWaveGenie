#ifndef __Positions__
#define __Positions__

#include <vector>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

class Positions
{
public:
    bool operator==(Positions& other);
    void insert(double x, double y, double z);
    typedef std::vector<double>::const_iterator ValueIterator;
    typedef boost::zip_iterator<boost::tuple<ValueIterator,ValueIterator,ValueIterator> > Iterator;
    size_t size();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator begin();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator end();
    void clear();
protected:
    std::vector<double> valuesX;
    std::vector<double> valuesY;
    std::vector<double> valuesZ;
};

#endif /* defined(__Positions__) */
