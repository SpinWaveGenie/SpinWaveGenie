#ifndef __GaussKronrod__
#define __GaussKronrod__

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <deque>
#include <memory>

class GaussKronrod
{
  //! Numerical Integration using the Gauss Kronrod quadrature.
  /*!
   This class vector quantities multiple dimensions.
   This class integrates a vector-valued integrand over multiple dimensions.

   The number of dimension is determined by the size of the vectors passed
   to setInterval and must match the number of evaluation points passed
   to m_integrand.

   Performance may not scale well to large numbers of dimensions.

   Method described on Wikipedia:
   http://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula
   */
public:
  GaussKronrod();
  GaussKronrod(const GaussKronrod &other);
  GaussKronrod &operator=(const GaussKronrod &other);
  GaussKronrod(GaussKronrod &&other);
  GaussKronrod &operator=(GaussKronrod &&other);
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
  //! \param maximumDivisions default is 1000
  void setMaximumDivisions(int maximumDivisions);
  //! Performs the integration
  //! \returns result of the integration.
  std::vector<double> integrate();
  ~GaussKronrod();

private:
  void setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints);
  class GKImpl;
  std::unique_ptr<GKImpl> m_p;
};

/*class GaussKronrod
{
public:
    GaussKronrod() : m_epsilon(1.0e-15),m_maxRecursionDepth(1.0e2) {};
    void setFunction(const std::function<
std::vector<double>(std::deque<double>& evaluationPoints)> & integrand)
    {
        m_integrand = integrand;
    };
    void setInterval(const std::vector<double>& lowerBounds, const
std::vector<double>& upperBounds);
    void setPrecision(double epsilon)
    {
        m_epsilon = epsilon;
    };
    void setMaximumRecursionDepth(long maxRecursionDepth)
    {
        m_maxRecursionDepth = maxRecursionDepth;
    };
    std::vector<double> integrate();
private:
    std::function< std::vector<double>(std::deque<double>& evaluationPoints)>
m_integrand;
    double m_lowerBound,m_upperBound,m_epsilon;
    std::vector<double> m_lowerBoundsInnerDimensions,
m_upperBoundsInnerDimensions;
    std::deque<double> m_evaluationPointsOuterDimensions;
    int m_maxRecursionDepth;
public:
    void setAdditionalEvaluationPoints(const std::deque<double>&
evaluationPoints)
    {
        m_evaluationPointsOuterDimensions = evaluationPoints;
    };
};
*/

#endif /* defined(__GaussKronrod__) */
