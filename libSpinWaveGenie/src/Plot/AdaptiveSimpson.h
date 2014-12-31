#ifndef __AdaptiveSimpsons__
#define __AdaptiveSimpsons__

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <deque>
#include <memory>

class AdaptiveSimpson
{
public:
    AdaptiveSimpson();
    AdaptiveSimpson(const AdaptiveSimpson& other);
    AdaptiveSimpson& operator= (const AdaptiveSimpson& other);
    AdaptiveSimpson(AdaptiveSimpson&& other);
    AdaptiveSimpson& operator= (AdaptiveSimpson&& other);
    void setFunction(const std::function< std::vector<double>(std::deque<double>& evaluationPoints)> & integrand);
    void setInterval(const std::vector<double>& lowerBounds, const std::vector<double>& upperBounds);
    void setPrecision(double epsilon);
    void setMaximumRecursionDepth(int maxRecursionDepth);
    std::vector<double> integrate();
    void setAdditionalEvaluationPoints(const std::deque<double>& evaluationPoints);
    ~AdaptiveSimpson();
private:
    class SimpsonImpl;
    std::unique_ptr<SimpsonImpl> m_p;
};

/*class AdaptiveSimpson
{
public:
    AdaptiveSimpson() : m_epsilon(1.0e-15),m_maxRecursionDepth(1.0e2) {};
    void setFunction(const std::function< std::vector<double>(std::deque<double>& evaluationPoints)> & integrand)
    {
        m_integrand = integrand;
    };
    void setInterval(const std::vector<double>& lowerBounds, const std::vector<double>& upperBounds);
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
    std::function< std::vector<double>(std::deque<double>& evaluationPoints)> m_integrand;
    double m_lowerBound,m_upperBound,m_epsilon;
    std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
    std::deque<double> m_evaluationPointsOuterDimensions;
    int m_maxRecursionDepth;
public:
    void setAdditionalEvaluationPoints(const std::deque<double>& evaluationPoints)
    {
        m_evaluationPointsOuterDimensions = evaluationPoints;
    };
};
*/

#endif /* defined(__AdaptiveSimpsons__) */

