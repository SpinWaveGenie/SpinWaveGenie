//
//  HKLDirection.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/22/14.
//
//

#ifndef __spin_wave_genie__HKLDirection__
#define __spin_wave_genie__HKLDirection__

#include <iostream>
#include <vector>

namespace SpinWaveGenie
{

//! Stores the integration direction and distance.
/*!
 This struct stores the integration direction and distance. v0,v1, and v2 fractions of H,K,and L
 and delta describes the integration distance.
 */
struct Axis
{
  //! Fraction of first reciprocal lattice vector in integration direction.
  double v0;
  //! Fraction of second reciprocal lattice vector in integration direction.
  double v1;
  //! Fraction of third reciprocal lattice vector in integration direction.
  double v2;
  //! Integration distance as a fraction of the vector shown above. Total distance is 2*delta
  double delta;
};

//! Stores the integration directions and distances used by IntegrateAxes.
/*!
 This container stores the integration directions and distances used by the IntegrateAxes class.
 It allows the user to integrate the calculated \f$ S\left(Q,\omega\right) \f$ over several unseen directions.
 Axes DO NOT need to be orthogonal.
 */

class HKLDirections
{
public:
  //! Add integration direction to container.
  //! \param direction Reciprocal lattice vector to integrate along. Acceptable inputs are 0,1, or 2.
  //! \param delta Integration distance as a fraction of the reciprocal lattice vector. Total distance is 2*delta.
  void addDirection(int direction, double delta);
  //! Add integration direction to container.
  //! \param v0 Fraction of first reciprocal lattice vector in integration direction.
  //! \param v1 Fraction of second reciprocal lattice vector in integration direction.
  //! \param v2 Fraction of third reciprocal lattice vector in integration direction.
  //! \param delta Integration distance as a fraction of the vector shown above. Total distance is 2*delta
  void addDirection(double v0, double v1, double v2, double delta);
  //! \return Number of integration directions in container.
  size_t size() const;
  //! Allow access to container via the subscript operator.
  //! \param position of Axis in container.
  //! \return reference to Axis object.
  const Axis &operator[](std::size_t position);
  typedef std::vector<Axis>::iterator Iterator;
  typedef std::vector<Axis>::const_iterator ConstIterator;
  //! \return Iterator pointing to the first Axis element
  Iterator begin();
  //! \return Iterator pointing to the end of the container.
  Iterator end();
  //! \return ConstIterator pointing to the first Axis element
  ConstIterator cbegin();
  //! \return ConstIterator pointing to the end of the container.
  ConstIterator cend();

private:
  std::vector<Axis> integrationDirections;
};

inline const Axis &HKLDirections::operator[](std::size_t position) { return integrationDirections[position]; }

inline size_t HKLDirections::size() const { return integrationDirections.size(); }

inline HKLDirections::Iterator HKLDirections::begin() { return integrationDirections.begin(); }

inline HKLDirections::Iterator HKLDirections::end() { return integrationDirections.end(); }

inline HKLDirections::ConstIterator HKLDirections::cbegin() { return integrationDirections.cbegin(); }

inline HKLDirections::ConstIterator HKLDirections::cend() { return integrationDirections.cend(); }
}
#endif /* defined(__spin_wave_genie__HKLDirection__) */
