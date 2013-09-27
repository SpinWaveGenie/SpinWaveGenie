#include "Neighbors.h"
#include "Cell.h"
#include "AtomIterator.h"
#include <boost/make_shared.hpp>

using namespace Eigen;
using namespace std;

Neighbors::Neighbors(boost::shared_ptr<Cell>& cell, string& sl1, string& sl2, double min, double max)
{
    neighborList = cell->getNeighbors(sl1, sl2, min, max);
}

AtomIterator Neighbors::begin()
{
    return AtomIterator(neighborList->begin());
}

AtomIterator Neighbors::end()
{
    return AtomIterator(neighborList->end());
}



