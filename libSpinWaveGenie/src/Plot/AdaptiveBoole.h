#ifndef __AdaptiveBoole__
#define __AdaptiveBoole__

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <deque>
#include <memory>

class AdaptiveBoole
//! Numerical Integration using the Adaptive Simpson's method.
/*!
 This class vector quantities multiple dimensions.
 This class integrates a vector-valued integrand over multiple dimensions.

 The number of dimension is determined by the size of the vectors passed
 to setInterval and must match the number of evaluation points passed
 to m_integrand.

 Performance may not scale well to large numbers of dimensions.

 Method described on Wikipedia:
 http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
 */
{
public:
  AdaptiveBoole();
  AdaptiveBoole(const AdaptiveBoole &other);
  AdaptiveBoole &operator=(const AdaptiveBoole &other);
  AdaptiveBoole(AdaptiveBoole &&other);
  AdaptiveBoole &operator=(AdaptiveBoole &&other);
  //! set function calculating the integrand.
  //! \param integrand function object must be of this type.
  void setFunction(const std::function<std::vector<double>(std::deque<double> &evaluationPoints)> &integrand);
  //! set interval in which to integrate over.
  //! \param lowerBounds array containing the lower bounds of each integral
  //! \param upperBounds array containing the upper bounds of each integral
  void setInterval(const std::vector<double> &lowerBounds, const std::vector<double> &upperBounds);
  //! sets the minimum estimated error.
  //! \param epsilon default is 1.0e-5.
  void setPrecision(double epsilon);
  //! sets the maximum number of times the algorithm with subdivide before returning.
  //! \param maximumDivisions default is 1000.
  void setMaximumDivisions(int maximumDivisions);
  //! Performs the integration
  //! \returns result of the integration.
  std::vector<double> integrate();
  ~AdaptiveBoole();

private:
  void setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints);
  class BooleImpl;
  std::unique_ptr<BooleImpl> m_p;
};

#endif /* defined(__AdaptiveBoole__) */
