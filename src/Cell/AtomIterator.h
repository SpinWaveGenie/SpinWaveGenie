#ifndef __AtomIterator_H__
#define __AtomIterator_H__

#include <vector>
#include <boost/iterator/iterator_facade.hpp>

/*!
Iterates over a number of atomic positions, such 
as those provided by the Sublattice or Neighbor class.
Example:
\code{.cpp}
Neighbors neighborList(cell,sl_rp,sl_sp,min,max);
AtomIterator nbrBegin = neighborList.begin();
AtomIterator nbrEnd = neighborList.end();
for(AtomIterator nbr=nbrBegin;nbr!=nbrEnd;++nbr)
{
    ...
}
\endcode
*/
class AtomIterator
: public boost::iterator_facade<
AtomIterator
, std::vector<double>
, boost::forward_traversal_tag
>
{
private:
    friend class boost::iterator_core_access;
    
public:
    explicit AtomIterator(std::vector<std::vector<double> >::iterator _it)
    : it(_it)
    {}
    
    AtomIterator(AtomIterator const& other)
    : it(other.it)
    {}
    
private:
    std::vector<std::vector<double> >::iterator it;
    
    void increment()
    {
        ++it;
    }
    
    bool equal(AtomIterator const& other) const
    {
        return (it == other.it);
    }
    
    std::vector<double>& dereference() const
    {
        return *it;
    }
    
};

#endif // __AtomIterator_H__