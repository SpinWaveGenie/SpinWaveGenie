#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <deque>
#include <functional>

class AdaptiveSimpson
{
    //! Numerical Integration using the Adaptive Simpson's method.
    /*!
     This class integrates a vector-valued integrand over multiple dimensions.
     
     The number of dimension is determined by the size of the vectors passed
     to setInterval and must match the number of evaluation points passed 
     to m_integrand.
     
     Performance may not scale well to large numbers of dimensions.
     
     Method described on Wikipedia:
     http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
     */
public:
    //! Default constructor set default value of epsilon and the recusion depth.
    AdaptiveSimpson() : m_epsilon(1.0e-5),m_maxRecursionDepth(10) {};
    //! set function calculating the integrand.
    //! \param integrand function object must be of this type.
    void setFunction(std::function< std::vector<double>(std::deque<double>& evaluationPoints)> const& integrand)
    {
        m_integrand = integrand;
    };
    //! set interval in which to integrate over.
    //! \param lowerBounds array containing the lower bounds of each integral
    //! \param upperBounds array containign the upper bounds of each integral
    void setInterval(std::vector<double>& lowerBounds, std::vector<double>& upperBounds);
    //! sets the minimum estimated error.
    //! \param epsilon default is 1.0e-5.
    void setPrecision(double epsilon)
    {
        m_epsilon = epsilon;
    };
    //! sets the maximum number of times the algorithm with subdivide before returning.
    //! \param maxRecursionDepth
    void setMaximumRecursionDepth(int maxRecursionDepth)
    {
        m_maxRecursionDepth = maxRecursionDepth;
    };
    //! Performs the integration
    //! \returns result of the integration.
    std::vector<double> integrate();
private:
    std::function< std::vector<double>(std::deque<double>& evaluationPoints)> m_integrand;
    double m_lowerBound,m_upperBound,m_epsilon;
    std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
    std::deque<double> m_evaluationPointsOuterDimensions;
    int m_maxRecursionDepth;
    std::vector<double> adaptiveSimpsons(double lowerBound, double upperBound, double epsilon,
                                         std::vector<double> S, std::vector<double>& fa,
                                         std::vector<double>& fb, std::vector<double>& fc,
                                         int recursionLevel);
    void setAdditionalEvaluationPoints(std::deque<double> evaluationPoints)
    {
        m_evaluationPointsOuterDimensions = evaluationPoints;
    };
};


