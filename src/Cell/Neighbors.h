#ifndef __Neighbors_H__
#define __Neighbors_H__

#define _USE_MATH_DEFINES
#include <iostream>
#include"Containers/Matrices.h"
#include "Containers/UniqueThreeVectors.h"

class Sublattice;
class Cell;

//! Finds neighbors of two sublattices between distances min and max.

class Neighbors
{
public:
    Neighbors();
    //! Finds neighbors of two sublattices between distances min and max.
    //! \param cell pointer to unit cell
    //! \param sl1 pointer to first sublattice
    //! \param sl2 pointer to second sublattice
    //! \param min minimum distance considered (Angstroms)
    //! \param max maximum distance considered (Angstroms)
    void findNeighbors(Cell& cell,std::string& sl1, std::string& sl2, double min, double max);
    //! Get the number of neighbors.
    double getNumberNeighbors();
    //! Get Gamma (descibed in paper J Phys. Condens. Matter 21 216001 (2009)
    //! \param K K vector used in spin wave calculation.
    std::complex<double> getGamma(Vector3 K);
    typedef UniqueThreeVectors<double>::Iterator Iterator;
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator begin();
    //! \return Returns an iterator pointing to the first element of the neighbor list
    Iterator end();
private:
    UniqueThreeVectors<double> neighborList;
    double numberNeighbors;
};
#endif // __Neighbors_H__ 
