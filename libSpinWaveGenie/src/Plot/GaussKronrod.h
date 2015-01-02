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
public:
  GaussKronrod();
  GaussKronrod(const GaussKronrod &other);
  GaussKronrod &operator=(const GaussKronrod &other);
  GaussKronrod(GaussKronrod &&other);
  GaussKronrod &operator=(GaussKronrod &&other);
  void setFunction(const std::function<std::vector<double>(std::deque<double> &evaluationPoints)> &integrand);
  void setInterval(const std::vector<double> &lowerBounds, const std::vector<double> &upperBounds);
  void setPrecision(double epsilon);
  void setMaximumRecursionDepth(int maxRecursionDepth);
  std::vector<double> integrate();
  void setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints);
  ~GaussKronrod();
private:
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
