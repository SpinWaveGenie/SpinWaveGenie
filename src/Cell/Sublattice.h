#ifndef __Sublattice_H__
#define __Sublattice_H__

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/iterator_facade.hpp>

class AtomIterator;

//! Describes a sublattice in the unit cell.
/*!
 The Sublattice class stores the orientation of the magnetic moment and
 the atomic positions of a particular sublattice. In addition, it calculates rotation
 and inverse rotation matrices.
 */
class Sublattice
{
public:
    //! Initializer
    Sublattice() {};
    //! Destructor
    ~Sublattice() {};
    //! use AtomIterator to iterate over atoms in sublattice
    //! set name to describe sublattice
    //! \param name_input a std::string argument
    void setName(std::string nameInput );
    //! returns name of a given sublattice
    //! \return name of sublattice
    std::string getName();
    //! set moment in spherical coordinates r,theta,phi
    /*! \param spin_input magnitude of spin moment
     \param theta_input angle 0 <= theta <= pi
     \param phi_input angle 0 <= phi <= 2*pi
     */
    void setMoment(double spinInput, double thetaInput , double phiInput);
    //! returns spherical coordinates (r,theta,phi) as a std:vector<double>
    //! \return coordinates [r,theta,phi]
    std::vector<double>* getMoment();
    //! returns rotation matrix as an Eigen::Matrix3d object
    //! \return rotation matrix
    Eigen::Matrix3d* getRotationMatrix();
    //! returns inverse rotation matrix as an Eigen::Matrix3d object
    //! \return inverse rotation matrix
    Eigen::Matrix3d* getInverseMatrix();
    //! add atom to the sublattice
    //! \param x x component of atomic position in Angstroms
    //! \param y y component of atomic position in Angstroms
    //! \param z z component of atomic position in Angstroms
    void addAtom(double x, double y, double z);
    //! returns iterator to first atomic position;
    AtomIterator begin();
    //! returns iterator to final atomic position;
    AtomIterator end();
private:
    std::string name; //!< name describing the sublattice
    double spin, //!< magnitude of spin
    theta, //!< angle theta describing orientation of spin
    phi; //!< angle phi describing orientation of spin
    std::vector<double> angles; //!< spin, theta, phi
    Eigen::Matrix3d rotationMatrix, //!< rotation matrix describing moment along z
    inverseMatrix; //!< inverse for rotation matrix
    std::vector<std::vector<double> > position; //!< std::vector storing atomic positions.
};
#endif // __Sublattice_H__