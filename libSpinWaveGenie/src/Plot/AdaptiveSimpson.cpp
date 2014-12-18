#include "AdaptiveSimpson.h"
#include <iterator>
#include <algorithm>

void AdaptiveSimpson::setInterval(std::vector<double>& lowerBounds, std::vector<double>& upperBounds)
{
    assert(lowerBounds.size() == upperBounds.size());
    m_lowerBound = lowerBounds.back();
    m_upperBound = upperBounds.back();
    if (lowerBounds.size()>1)
    {
        m_lowerBoundsInnerDimensions = lowerBounds;
        m_lowerBoundsInnerDimensions.pop_back();
        
        m_upperBoundsInnerDimensions = upperBounds;
        m_upperBoundsInnerDimensions.pop_back();
    }
};

std::vector<double> AdaptiveSimpson::integrate()
{
    double c = (m_lowerBound + m_upperBound)/2.0;
    std::vector<double> fa,fb,fc;
    if (m_lowerBoundsInnerDimensions.size() > 0)
    {
        AdaptiveSimpson test;
        test.setFunction(m_integrand);
        test.setInterval(m_lowerBoundsInnerDimensions,m_upperBoundsInnerDimensions);
        test.setMaximumRecursionDepth(m_maxRecursionDepth);
        test.setPrecision(m_epsilon);
        m_evaluationPointsOuterDimensions.push_front(m_lowerBound);
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        fa = test.integrate();
        m_evaluationPointsOuterDimensions[0] = m_upperBound;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        fb = test.integrate();
        m_evaluationPointsOuterDimensions[0] = c;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        fc = test.integrate();
    }
    else
    {
        m_evaluationPointsOuterDimensions.push_front(m_lowerBound);
        fa = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = m_upperBound;
        fb = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = c;
        fc = m_integrand(m_evaluationPointsOuterDimensions);
    }
    std::vector<double> S;
    std::size_t size = fa.size();
    assert(fb.size() == size);
    assert(fc.size() == size);
    S.reserve(size);
    double prefactor = (m_upperBound - m_lowerBound)/6.0;
    for (std::size_t i = 0; i<size; i++)
    {
        S.push_back(prefactor*(fa[i] + 4.0*fc[i] + fb[i]));
    }
    return this->adaptiveSimpsons(m_lowerBound, m_upperBound, m_epsilon, S, fa, fb, fc,0);
};


double calculateResult(double S,double S2)
{
    return S2 + (S2 - S)/15.0;
}

//
// Recursive auxiliary function for adaptiveSimpsons() function below
//
std::vector<double> AdaptiveSimpson::adaptiveSimpsons(double lowerBound, double upperBound, double epsilon,
                                                       std::vector<double> S, std::vector<double>& fa, std::vector<double>& fb, std::vector<double>& fc, int recursionLevel) {
    double c = (lowerBound + upperBound)/2, h = upperBound - lowerBound;
    double d = (lowerBound + c)/2, e = (c + upperBound)/2;
    std::vector<double> fd,fe;
    if (m_lowerBoundsInnerDimensions.size() > 0)
    {
        AdaptiveSimpson test;
        test.setFunction(m_integrand);
        test.setInterval(m_lowerBoundsInnerDimensions,m_upperBoundsInnerDimensions);
        test.setMaximumRecursionDepth(m_maxRecursionDepth);
        test.setPrecision(m_epsilon);
        m_evaluationPointsOuterDimensions[0] = d;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        fd = test.integrate();
        m_evaluationPointsOuterDimensions[0] = e;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        fe = test.integrate();
    }
    else
    {
        m_evaluationPointsOuterDimensions[0] = d;
        fd = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = e;
        fe = m_integrand(m_evaluationPointsOuterDimensions);
        
    }
    std::size_t size = fd.size();
    assert(fe.size() == size);
    std::vector<double> Sleft,Sright,S2;
    Sleft.reserve(size);Sright.reserve(size); S2.reserve(size);
    bool done = true;
    double prefactor = h/12.0;
    double eps_comparison = 15.0*epsilon;
    for (std::size_t i = 0; i< size; i++)
    {
        Sleft.push_back(prefactor*(fa[i] + 4.0*fd[i] + fc[i]));
        Sright.push_back(prefactor*(fc[i] + 4.0*fe[i] + fb[i]));
        S2.push_back(Sleft[i] + Sright[i]);
        if (done && std::abs(S2[i] - S[i]) <= eps_comparison)
            done = true;
        else
            done = false;
    }
    if (recursionLevel > m_maxRecursionDepth)
    {
        done = true;
    }
    
    std::vector<double> result;
    result.reserve(size);
    if (done)
    {
        std::insert_iterator<std::vector<double> > insertResult(result,result.begin());
        std::transform (S.begin(), S.end(), S2.begin(), insertResult, calculateResult);
    }
    else
    {
        std::vector<double> left = adaptiveSimpsons(lowerBound, c, epsilon/2.0, Sleft,  fa, fc, fd, recursionLevel+1);
        std::vector<double> right = adaptiveSimpsons(c, upperBound, epsilon/2.0, Sright, fc, fb, fe, recursionLevel+1);
        
        std::insert_iterator<std::vector<double> > insertResult(result,result.begin());
        std::transform (left.begin(), left.end(), right.begin(), insertResult, std::plus<double>());
    }
    return result;
}
