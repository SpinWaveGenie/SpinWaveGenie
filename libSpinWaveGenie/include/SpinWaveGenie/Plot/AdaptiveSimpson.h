#ifndef __AdaptiveSimpson__
#define __AdaptiveSimpson__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Memory.h"
#include <deque>
#include <functional>
#include <vector>

class SPINWAVEGENIE_EXPORT AdaptiveSimpson
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
  AdaptiveSimpson();
  AdaptiveSimpson(const AdaptiveSimpson &other);
  AdaptiveSimpson &operator=(const AdaptiveSimpson &other);
  AdaptiveSimpson(AdaptiveSimpson &&other);
  AdaptiveSimpson &operator=(AdaptiveSimpson &&other);
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
  void setMaximumDivisions(std::size_t maximumDivisions);
  //! Performs the integration
  //! \returns result of the integration.
  std::vector<double> integrate();
  ~AdaptiveSimpson();

private:
  void setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints);
  class SimpsonImpl;
  std::unique_ptr<SimpsonImpl> m_p;
};

#endif /* defined(__AdaptiveSimpson__) */
