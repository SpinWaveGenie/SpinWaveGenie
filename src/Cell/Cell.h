//
#ifndef __Cell_H__
#define __Cell_H__

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Sublattice.h"
#include "Containers/Matrices.h"

//! Unit cell containing basis vectors and all sublattices.
/*!
The Cell class stores the basis vectors and all sublattices in the unit cell.
Atomic positions are inserted as a fraction of the basis vectors and converted to Angstroms.
*/
class Cell
{
public:
    //! Set basis vectors from parameters a, b, c, \f$\alpha,\: \beta,\: \gamma \f$.
    //! \param a     Distance a in Angstroms
    //! \param b     Distance b in Angstroms
    //! \param c     Distance c in Angstroms
    //! \param alpha Angle \f$ \alpha \f$ in radians
    //! \param beta  Angle \f$ \beta \f$ in radians
    //! \param gamma Angle \f$ \gamma \f$ in radians
    void setBasisVectors(double a,double b, double c, double alpha, double beta, double gamma);
    //! Set basis vectors as Eigen::Matrix3d object. Vectors are stored as rows.
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
    //! \param sl Sublattice object
    //! \return reciprocal lattice vectors as rows in an Eigen::Matrix3d object
    void addSublattice(Sublattice& sl);
    //! Returns sublattice "name"
    //! \param name used to describe sublattice
    //! \return sublattice
    Sublattice& getSublattice(std::string name);
    //! Add atom to sublattice name at position pos
    //! \param name Sublattice atom belongs to
    //! \param x x coordinate of atom in fraction of the basis vectors.
    //! \param y y coordinate of atom in fraction of the basis vectors.
    //! \param z z coordinate of of atom in fraction of the basis vectors.
    void addAtom(std::string name, double x, double y, double z);
    //! Returns the position where sublattice name is stored.
    //! \param name name of sublattice.
    const std::size_t getPosition(std::string name);
    //! Returns the number of sublattices in the cell
    //! \return number of sublattices
    const size_t size() const;
    typedef std::vector<Sublattice>::iterator Iterator;
    typedef std::vector<Sublattice>::const_iterator ConstIterator;
    //! \return Returns an Iterator pointing to the first Sublattice element
    Iterator begin();
    //! \return Returns an Iterator pointing to the end of the vector
    Iterator end();
    //! \return Returns a ConstIterator pointing to the first Sublattice element
    ConstIterator cbegin();
    //! \return Returns a ConstIterator pointing to the end of the vector
    ConstIterator cend();
private:
    Matrix3 basisVectors;
    Matrix3 reciprocalVectors;
    std::vector<Sublattice> sublatticeInfo;
};
#endif // __Cell_H__ 
