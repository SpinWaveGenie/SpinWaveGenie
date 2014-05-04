//
#ifndef __Cell_H__
#define __Cell_H__

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <memory> 
#include <boost/iterator/iterator_facade.hpp>
#include "Sublattice.h"
#include "Containers/Matrices.h"

//! Contains the basis vectors and all sublattices in the unit cell.
/*!
The Cell class stores the basis vectors and all sublattices in the unit cell.
Atomic positions are inserted as a fraction of the basis vectors and converted to Angstroms.
*/
class Cell
{
public:
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
    void setBasisVectors(double scale, Matrix3 basis);
    //! get basis vectors as an Eigen::Matrix3d object
    //! \return basis vectors as rows in an Eigen::Matrix3d object
    const Matrix3 getBasisVectors() const;
    //! get basis vectors as an Eigen::Matrix3d object
    //! \return basis vectors as rows in an Eigen::Matrix3d object
    const Matrix3 getReciprocalVectors() const;
    //! Add sublattice to cell
    //! \param name Unique name defining sublattice
    //! \param sl Sublattice object
    //! \return reciprocal lattice vectors as rows in an Eigen::Matrix3d object
    void addSublattice(Sublattice& sl);
    //! Returns sublattice "name"
    //! \param name used to describe sublattice
    //! \return sublattice
    Sublattice& getSublattice(std::string name) const;
    //! Add atom to sublattice name at position pos
    //! \param name Sublattice atom belongs to
    //! \param pos Position of atom in fraction of the basis vectors.
    const std::size_t getPosition(std::string name);
    void addAtom(std::string name, double x, double y, double z);
    //! Returns the number of sublattices in the cell
    //! \return number of sublattices
    const size_t size() const;
    typedef std::vector<Sublattice>::iterator Iterator;
    typedef std::vector<Sublattice>::const_iterator ConstIterator;
    //! \return Returns an iterator pointing to the first Sublattice element
    Iterator begin();
    ConstIterator cbegin();
    //! \return Returns an iterator pointing to the final Sublattice element
    Iterator end();
    ConstIterator cend();
    //~Cell() { std::cout << "Cell destructed" << std::endl; }
private:
    //! basis vectors
    Matrix3 basisVectors;
    //! reciprocal lattice vectors
    Matrix3 reciprocalVectors;
    //! vector containing all Sublattice objects;
    std::vector<Sublattice> sublatticeInfo;
};
#endif // __Cell_H__ 
