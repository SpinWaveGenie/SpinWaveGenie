#ifndef __Sublattice_H__
#define __Sublattice_H__

#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"
#include "SpinWaveGenie/Export.h"
#include <string>
#include <vector>

namespace SpinWaveGenie
{

//! Sublattice in the unit cell.
/*!
 The Sublattice class stores the orientation of the magnetic moment and
 the atomic positions of a particular sublattice. In addition, it calculates rotation
 and inverse rotation matrices.
 */
class SPINWAVEGENIE_EXPORT Sublattice
{
public:
  //! Default constructor
  Sublattice();
  //! set name to describe sublattice
  //! \param nameInput a std::string argument unique to each sublattice
  void setName(const std::string &nameInput);
  //! returns name of a given sublattice
  //! \return name of sublattice
  std::string getName() const;
  //! set type to describe magnetic form factor used in the calculation of intensities
  //! \param typeInput a std::string argument
  void setType(const std::string &typeInput);
  //! returns name of a given sublattice
  //! \return name of sublattice
  std::string getType() const;
  //! set moment in spherical coordinates r,theta,phi
  /*! \param spinInput magnitude of spin moment
   \param thetaInput angle \f$ 0 \leq \theta \leq  \pi \f$
   \param phiInput angle \f$ 0 \leq \phi \leq 2\pi \f$
   */
  void setMoment(double spinInput, double thetaInput, double phiInput);
  //! \return coordinate \f$ r  \f$ of \f$ \left( r,\theta,\phi \right) \f$, unitless
  double getMoment() const { return spin; }
  //! \return coordinates \f$\theta\f$ of \f$ \left( r,\theta,\phi \right) \f$, in radians
  double getTheta() const { return theta; }
  //! \return coordinates \f$ \phi \f$ of \f$ \left( r,\theta,\phi \right) \f$, in radians
  double getPhi() const { return phi; }
  //! returns rotation matrix as an Eigen::Matrix3d object
  //! \return rotation matrix
  const Matrix3 &getRotationMatrix() const { return rotationMatrix; }
  //! returns inverse rotation matrix as an Eigen::Matrix3d object
  //! \return inverse rotation matrix
  const Matrix3 &getInverseMatrix() const { return inverseMatrix; }
  //! add atom to the sublattice
  //! \param x x component of atomic position in Angstroms
  //! \param y y component of atomic position in Angstroms
  //! \param z z component of atomic position in Angstroms
  void addAtom(double x, double y, double z);
  //! returns the number of atoms in this sublattice
  //! \return number of atoms
  std::size_t size() const { return positions.size(); }
  typedef UniqueThreeVectors<double>::Iterator Iterator;
  typedef UniqueThreeVectors<double>::ConstIterator ConstIterator;
  //! returns an Iterator to the first atomic position
  Iterator begin() { return positions.begin(); }
  //! returns an Iterator to the end of the vector
  Iterator end() { return positions.end(); }
  //! returns an Iterator to the first atomic position
  ConstIterator begin() const { return positions.cbegin(); }
  //! returns an Iterator to the end of the vector
  ConstIterator end() const { return positions.cend(); }
  //! returns a ConstIterator to the first atomic position
  ConstIterator cbegin() const { return positions.cbegin(); }
  //! returns a ConstIterator to the end of the vector
  ConstIterator cend() const { return positions.cend(); }

private:
  std::string name, type;
  double spin, theta, phi;
  Matrix3 rotationMatrix, inverseMatrix;
  UniqueThreeVectors<double> positions;
};
}

#endif // __Sublattice_H__
