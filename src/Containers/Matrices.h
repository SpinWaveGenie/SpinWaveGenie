//
//  EigenStuff.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/18/13.
//
//
#ifndef __spin_wave_genie__Matrices__
#define __spin_wave_genie__Matrices__

#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "Eigen/Core"

using Eigen::Dynamic;
using Eigen::RowMajor;

typedef Eigen::Matrix <double, 3, 1> Vector3;
typedef Eigen::Matrix<double, Dynamic, 1> Vector;

typedef Eigen::Matrix <double, 3, 3> Matrix3;
typedef Eigen::Matrix<double, Dynamic, Dynamic> Matrix;
typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> MatrixRowMajor;

/*!
 Iterates over a number of atomic positions, such
 as those provided by the Sublattice or Neighbor class.
 Example:
 \code{.cpp}
 Neighbors neighborList(cell,sl_rp,sl_sp,min,max);
 for(Neighbors::Iterator nbr=neighborList.begin();nbr!=neighborList.end();++nbr)
 {
 ...
 }
 \endcode
 */

class Vector3Iterator
: public boost::iterator_facade<
Vector3Iterator
, Vector3
, boost::forward_traversal_tag
>
{
private:
    friend class boost::iterator_core_access;
    
public:
    explicit Vector3Iterator(std::vector<Vector3>::iterator _it)
    : it(_it)
    {}
    
    Vector3Iterator( Vector3Iterator const& other)
    : it(other.it)
    {}
    
private:
    std::vector<Vector3>::iterator it;
    
    void increment()
    {
        ++it;
    }
    
    bool equal(Vector3Iterator const& other) const
    {
        return (it == other.it);
    }
    
    Vector3& dereference() const
    {
        return *it;
    }
    
};


#endif /* defined(__spin_wave_genie__Matrices__) */
