#ifndef __Sublattice_H__
#define __Sublattice_H__

#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Containers/Matrices.h"
#include "Containers/UniqueThreeVectors.h"

//! Describes a sublattice in the unit cell.
/*!
 The Sublattice class stores the orientation of the magnetic moment and
 the atomic positions of a particular sublattice. In addition, it calculates rotation
 and inverse rotation matrices.
 */
class Sublattice
{
public:
    Sublattice();
    ~Sublattice() {};
    //! set name to describe sublattice
    //! \param name_input a std::string argument
    void setName(std::string nameInput );
    //! returns name of a given sublattice
    //! \return name of sublattice
    std::string getName() const;
    //! set type to describe magnetic form factor used in the calculation of intensities;
    //! \param typeInput a std::string argument
    void setType(std::string typeInput );
    //! returns name of a given sublattice
    //! \return name of sublattice
    std::string getType() const;
    //! set moment in spherical coordinates r,theta,phi
    /*! \param spin_input magnitude of spin moment
     \param theta_input angle 0 <= theta <= pi
     \param phi_input angle 0 <= phi <= 2*pi
     */
    void setMoment(double spinInput, double thetaInput , double phiInput);
    //! \return coordinate r of [r,theta,phi]
    double getMoment() const;
    //! \return coordinates theta of [r,theta,phi]
    double getTheta() const;
    //! \return coordinates phi of [r,theta,phi]
    double getPhi() const;
    //! returns rotation matrix as an Eigen::Matrix3d object
    //! \return rotation matrix
    const Matrix3& getRotationMatrix() const;
    //! returns inverse rotation matrix as an Eigen::Matrix3d object
    //! \return inverse rotation matrix
    const Matrix3& getInverseMatrix() const;
    //! add atom to the sublattice
    //! \param x x component of atomic position in Angstroms
    //! \param y y component of atomic position in Angstroms
    //! \param z z component of atomic position in Angstroms
    void addAtom(double x, double y, double z);
    typedef UniqueThreeVectors<double>::Iterator Iterator;
    typedef UniqueThreeVectors<double>::ConstIterator ConstIterator;
    //! returns iterator to first atomic position;
    Iterator begin();
    //! returns iterator to final atomic position;
    Iterator end();
    //! returns iterator to first atomic position;
    ConstIterator cbegin();
    //! returns iterator to final atomic position;
    ConstIterator cend();
private:
    std::string name; //!< name describing the sublattice
    std::string type; //!< type describing the magnetic form factor
    double spin, //!< magnitude of spin
    theta, //!< angle theta describing orientation of spin
    phi; //!< angle phi describing orientation of spin
    Matrix3 rotationMatrix, //!< rotation matrix describing moment along z
    inverseMatrix; //!< inverse for rotation matrix
    UniqueThreeVectors<double> positions; //!< std::vector storing atomic positions.
};
#endif // __Sublattice_H__
