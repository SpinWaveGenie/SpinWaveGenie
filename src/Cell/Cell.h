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
    //! Use CellIter to iterate over sublattices in the unit cell
    friend class CellIter;
    //! Use Neighbors to iterate over the relative position of neighboring atoms of a particular sublattice
    friend class Neighbors;
    //! Set basis vectors as Eigen::Matrix3d object. May be multiplied by double "scale"
    //! \param scale Double used to scale basis vectors
    //! \param basis Basis vectors as rows in an Eigen::Matrix3d object
    void set_basis_vectors(double a,double b, double c, double alpha, double beta, double gamma);
    void set_basis_vectors(double scale, Eigen::Matrix3d basis);
    //! get basis vectors as an Eigen::Matrix3d object
    //! \return basis vectors as rows in an Eigen::Matrix3d object
    Eigen::Matrix3d get_basis_vectors();
    Eigen::Matrix3d get_reciprocal_vectors();
    //! Add sublattice to cell
    //! \param name
    //! \param sl
    void add_sublattice(std::string name, boost::shared_ptr<Sublattice>& sl);
    //! Returns pointer to sublattice "name"
    //! \param name used to describe sublattice
    //! \return pointer to sublattice
    boost::shared_ptr<Sublattice> get_sublattice(std::string name);
    //! Add atom to sublattice name at position pos
    //! \param name Sublattice atom belongs to
    //! \param pos Position of atom in fraction of the basis vectors.
    void add_atom(std::string name, double x, double y, double z);
    //! Returns the number of sublattices in the cell
    // \return number of sublattices
    int size();
    SublatticeIterator begin();
    SublatticeIterator end();
private:
    struct FastCompare
    {
        boost::shared_ptr<Sublattice> sl1, sl2;
        double min,max;
        bool operator ==(const FastCompare &other) const;
        bool operator < ( const FastCompare &other) const;
    };
    //! Finds neighbors of two sublattices between distances min and max.
    //! Used exclusively by NeighborIter
    //! \param sublattice1_in pointer to first sublattice
    //! \param sublattice2_in pointer to second sublattice
    //! \param min minimum distance considered (Angstroms)
    //! \param max maximum distance considered (Angstroms)
    std::vector<std::vector<double> >* get_neighbors(boost::shared_ptr<Sublattice>& sublattice1_in, boost::shared_ptr<Sublattice>& sublattice2_in, double min, double max);
    // returns relative position of all neighboring atoms
    Eigen::Matrix3d basis_vectors,reciprocal_vectors;
    std::map<std::string,boost::shared_ptr<Sublattice> > sublattice_info;
    std::map<FastCompare,std::vector<std::vector<double> > > neighborCache;
};

#endif // __Cell_H__ 
