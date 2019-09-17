//
#ifndef __Cell_H__
#define __Cell_H__

#include "Eigen/Core"
#include "SpinWaveGenie/Containers/Sublattice.h"
#include "SpinWaveGenie/Export.h"
#include <string>
#include <vector>

namespace SpinWaveGenie
{
//! Unit cell containing basis vectors and all sublattices.
/*!
The Cell class stores the basis vectors and all sublattices in the unit cell.
Atomic positions are inserted as a fraction of the basis vectors and converted to Angstroms.
*/
class SPINWAVEGENIE_EXPORT Cell
{
public:
  //! Set basis vectors from parameters a, b, c, \f$\alpha,\: \beta,\: \gamma \f$.
  //! \param a     Distance a in Angstroms
  //! \param b     Distance b in Angstroms
  //! \param c     Distance c in Angstroms
  //! \param alpha Angle \f$ \alpha \f$ in degrees
  //! \param beta  Angle \f$ \beta \f$ in degrees
  //! \param gamma Angle \f$ \gamma \f$ in degrees
  void setBasisVectors(double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg);
  //! Set basis vectors as Eigen::Matrix3d object. Vectors are stored as rows.
  //! \param scale Double used to scale basis vectors
  //! \param basis Basis vectors as rows in an Eigen::Matrix3d object
  void setBasisVectors(double scale, const Eigen::Matrix3d &basis);
  //! get basis vectors as an Eigen::Matrix3d object
  //! \return basis vectors as rows in an Eigen::Matrix3d object
  const Eigen::Matrix3d &getBasisVectors() const;
  //! get basis vectors as an Eigen::Matrix3d object
  //! \return basis vectors as rows in an Eigen::Matrix3d object
  const Eigen::Matrix3d &getReciprocalVectors() const;
  //! Add sublattice to cell
  //! \param sl Sublattice object
  //! \return reciprocal lattice vectors as rows in an Eigen::Matrix3d object
  void addSublattice(const Sublattice &sl);
  //! Returns sublattice "name"
  //! \param name used to describe sublattice
  //! \return sublattice
  Sublattice &getSublattice(const std::string &name);
  const Sublattice &getSublattice(const std::string &name) const;
  const Sublattice &operator[](std::vector<Sublattice>::size_type position) const;
  //! Add atom to sublattice name at position pos
  //! \param name Sublattice atom belongs to
  //! \param x x coordinate of atom in fraction of the basis vectors.
  //! \param y y coordinate of atom in fraction of the basis vectors.
  //! \param z z coordinate of of atom in fraction of the basis vectors.
  void addAtom(const std::string &name, double x, double y, double z);
  //! Returns the position where sublattice name is stored.
  //! \param name name of sublattice.
  std::vector<Sublattice>::difference_type getPosition(const std::string &name) const;
  //! Returns the number of sublattices in the cell
  //! \return number of sublattices
  std::vector<Sublattice>::size_type size() const { return sublatticeInfo.size(); }
  using Iterator = std::vector<Sublattice>::iterator;
  using ConstIterator = std::vector<Sublattice>::const_iterator;
  //! \return Returns an Iterator pointing to the first Sublattice element
  Iterator begin() { return sublatticeInfo.begin(); }
  //! \return Returns an Iterator pointing to the end of the vector
  Iterator end() { return sublatticeInfo.end(); }
  //! \return Returns a ConstIterator pointing to the first Sublattice element
  ConstIterator begin() const { return sublatticeInfo.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the vector
  ConstIterator end() const { return sublatticeInfo.cend(); }
  //! \return Returns a ConstIterator pointing to the first Sublattice element
  ConstIterator cbegin() const { return sublatticeInfo.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the vector
  ConstIterator cend() const { return sublatticeInfo.cend(); }

private:
  Eigen::Matrix3d basisVectors;
  Eigen::Matrix3d reciprocalVectors;
  std::vector<Sublattice> sublatticeInfo;
};
} // namespace SpinWaveGenie
#endif // __Cell_H__
