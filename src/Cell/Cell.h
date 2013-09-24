//
//  Cell.h
//  Spin Wave Genie
//
//  Created by Hahn, Steven E. on 2/6/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//

#ifndef __Cell_H__
#define __Cell_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "Sublattice.h"

//! SublatticeIterator class
/*!
 The SublatticeIterator class iterates over sublattices stored in the Cell class. Example:
 \code{.cpp}
 int r;
 int M = 0;
 for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
 {
     if ( sl_r == (*sl)->getName())
     r = M;
     M++;
 }
 \endcode
 */
class SublatticeIterator
: public boost::iterator_facade<
SublatticeIterator
, boost::shared_ptr<Sublattice>
, boost::forward_traversal_tag
>
{
private:
    friend class boost::iterator_core_access;
    
public:
    explicit SublatticeIterator(std::map<std::string,boost::shared_ptr<Sublattice> >::iterator _it)
    : it(_it)
    {}
    
    SublatticeIterator(SublatticeIterator const& other)
    : it(other.it)
    {}
    
private:
    std::map<std::string,boost::shared_ptr<Sublattice> >::iterator it;
    
    void increment()
    {
        ++it;
    }
    
    bool equal(SublatticeIterator const& other) const
    {
        return (it == other.it);
    }
    
    boost::shared_ptr<Sublattice>& dereference() const
    {
        return it->second;
    }
    
};

//! Cell class
/*!
The Cell class stores the basis vectors and a pointer to all sublattices in the unit cell. 
Atomic positions inserted as a fraction of the basis vectors and converted to Angstroms.
*/
class Cell
{
public:
    //! Use Neighbors class to iterate over the relative position of neighboring atoms of a particular sublattice
    friend class Neighbors;
    //! Set basis vectors from variable a,b,c,alpha,beta,gamma.
    //! \param a     Distance a in Angstroms
    //! \param b     Distance b in Angstroms
    //! \param c     Distance c in Angstroms
    //! \param alpha Angle alpha in degrees
    //! \param beta  Angle beta in degrees
    //! \param gamma Angle gamma in degrees
    void setBasisVectors(double a,double b, double c, double alpha, double beta, double gamma);
    //! Set basis vectors as Eigen::Matrix3d object. May be multiplied by double "scale"
    //! \param scale Double used to scale basis vectors
    //! \param basis Basis vectors as rows in an Eigen::Matrix3d object
    void setBasisVectors(double scale, Eigen::Matrix3d basis);
    //! get basis vectors as an Eigen::Matrix3d object
    //! \return basis vectors as rows in an Eigen::Matrix3d object
    Eigen::Matrix3d getBasisVectors();
    //! get basis vectors as an Eigen::Matrix3d object
    //! \return basis vectors as rows in an Eigen::Matrix3d object
    Eigen::Matrix3d getReciprocalVectors();
    //! Add sublattice to cell
    //! \param name Unique name defining sublattice
    //! \param sl Pointer to Sublattice object
    //! \return reciprocal lattice vectors as rows in an Eigen::Matrix3d object
    void addSublattice(std::string name, boost::shared_ptr<Sublattice>& sl);
    //! Returns pointer to sublattice "name"
    //! \param name used to describe sublattice
    //! \return pointer to sublattice
    boost::shared_ptr<Sublattice> getSublattice(std::string name);
    //! Add atom to sublattice name at position pos
    //! \param name Sublattice atom belongs to
    //! \param pos Position of atom in fraction of the basis vectors.
    void addAtom(std::string name, double x, double y, double z);
    //! Returns the number of sublattices in the cell
    //! \return number of sublattices
    int size();
    //! \return Returns an iterator pointing to the first Sublattice element
    SublatticeIterator begin();
    //! \return Returns an iterator pointing to the final Sublattice element
    SublatticeIterator end();
private:
    //! Finds neighbors of two sublattices between distances min and max.
    //! Used exclusively by NeighborIter
    //! \param sublattice1_in pointer to first sublattice
    //! \param sublattice2_in pointer to second sublattice
    //! \param min minimum distance considered (Angstroms)
    //! \param max maximum distance considered (Angstroms)
    //! \return relative position of all neighboring atoms
    std::vector<std::vector<double> >* get_neighbors(boost::shared_ptr<Sublattice>& sublattice1_in, boost::shared_ptr<Sublattice>& sublattice2_in, double min, double max);
    //! Finding the relative position of neighboring atoms is a relatively slow process,
    //! so we want to store the results so they can be reused. Therefore, we need a struct that quickly identifies when the input parameters are identical.
    struct FastCompare
    {
        boost::shared_ptr<Sublattice> sl1, sl2;
        double min,max;
        bool operator ==(const FastCompare &other) const;
        bool operator < ( const FastCompare &other) const;
    };
    //! stores the results of getNeighbors in a map so they can be reused.
    std::map<FastCompare,std::vector<std::vector<double> > > neighborCache;
    //! basis vectors
    Eigen::Matrix3d basisVectors;
    //! reciprocal lattice vectors
    Eigen::Matrix3d reciprocalVectors;
    //! Pointers to all Sublattice objects;
    std::map<std::string,boost::shared_ptr<Sublattice> > sublatticeInfo;
};

#endif // __Cell_H__ 
