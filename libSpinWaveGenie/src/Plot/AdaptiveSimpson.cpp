#include "AdaptiveSimpson.h"

#include <queue>
#include <algorithm>

struct helper
{
    std::vector<double> fa,fb,fc,fd,fe;
    std::vector<double> S,Sleft,Sright;
    double lowerlimit;
    double upperlimit;
    double epsilon;
    double error = 0.0;
    bool operator<(const helper& rhs) const
    {
        return this->error/this->epsilon < rhs.error/rhs.epsilon;
    }
};

class AdaptiveSimpson::SimpsonImpl
{
public:
    std::vector<double> sumPieces(std::priority_queue<helper>& pieces);
    helper createElement(const helper& mostError,bool first);
    std::pair<helper,helper> splitElement(const helper& mostError);
    std::vector<double> integrate();
    std::function< std::vector<double>(std::deque<double>& evaluationPoints)> m_integrand;
    double m_lowerBound,m_upperBound,m_epsilon = 1.0e-3;
    std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
    std::deque<double> m_evaluationPointsOuterDimensions;
    int m_maxRecursionDepth = 1000;
    std::unique_ptr<SimpsonImpl> clone();
};

std::unique_ptr<AdaptiveSimpson::SimpsonImpl> AdaptiveSimpson::SimpsonImpl::clone()
{
    return std::unique_ptr<AdaptiveSimpson::SimpsonImpl>(new SimpsonImpl(*this));
}

helper AdaptiveSimpson::SimpsonImpl::createElement(const helper& mostError, bool first)
{
    helper element;
    double a,b;
    if (first)
    {
        a = mostError.lowerlimit;
        b = (mostError.lowerlimit+mostError.upperlimit)/2.0;
    }
    else
    {
        a = (mostError.lowerlimit+mostError.upperlimit)/2.0;
        b = mostError.upperlimit;
    }
    double c = (a+b)/2.0;
    double d = (a+c)/2.0;
    double e = (c+b)/2.0;
    if (first)
    {
        element.fa = std::move(mostError.fa);
        element.fb = mostError.fc;
        element.fc = std::move(mostError.fd);
        element.S = std::move(mostError.Sleft);
    }
    else
    {
        element.fa = std::move(mostError.fc);
        element.fb = std::move(mostError.fb);
        element.fc = std::move(mostError.fe);
        element.S = std::move(mostError.Sright);
    }
    
    if (m_lowerBoundsInnerDimensions.size() > 0)
    {
        AdaptiveSimpson test;
        test.setFunction(m_integrand);
        test.setInterval(m_lowerBoundsInnerDimensions,m_upperBoundsInnerDimensions);
        test.setMaximumRecursionDepth(m_maxRecursionDepth);
        test.setPrecision(m_epsilon);
        m_evaluationPointsOuterDimensions[0] = d;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        element.fd = test.integrate();
        m_evaluationPointsOuterDimensions[0] = e;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        element.fe = test.integrate();
    }
    else
    {
        m_evaluationPointsOuterDimensions[0] = d;
        element.fd = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = e;
        element.fe = m_integrand(m_evaluationPointsOuterDimensions);
    }
    
    std::size_t size = element.fa.size();

    double prefactor = (b-a)/12.0;
    element.Sleft.reserve(size);
    element.Sright.reserve(size);
    for (std::size_t i = 0; i<size; i++)
    {
        element.Sleft.emplace_back(prefactor*(element.fa[i] + 4.0*element.fd[i] + element.fc[i]));
        element.Sright.emplace_back(prefactor*(element.fc[i] + 4.0*element.fe[i] + element.fb[i]));
        element.error = std::max(element.error,std::abs(15.0*(element.Sleft[i]+element.Sright[i]-element.S[i])));
    }
    element.lowerlimit = a;
    element.upperlimit = b;
    element.epsilon = std::max(mostError.epsilon/sqrt(2.0),std::numeric_limits<double>::epsilon());
    return element;
}

std::pair<helper,helper> AdaptiveSimpson::SimpsonImpl::splitElement(const helper& mostError)
{
    helper element1 = this->createElement(mostError,true);
    helper element2 = this->createElement(mostError,false);
    return std::tie(element1,element2);
}

std::vector<double> AdaptiveSimpson::SimpsonImpl::sumPieces(std::priority_queue<helper>& pieces)
{
    std::size_t size = pieces.top().Sleft.size();
    std::vector<double> sum(size);
    while (pieces.size() > 0)
    {
        auto element = pieces.top();
        for (std::size_t i = 0; i < size; i++)
        {
            double S2 = element.Sleft[i] + element.Sright[i];
            sum[i] +=  S2 + (S2 - element.S[i])/15.0;
        }
        pieces.pop();
    }
    return sum;
}

std::vector<double> AdaptiveSimpson::SimpsonImpl::integrate()
{
    helper first;
    
    double a = m_lowerBound;
    double b = m_upperBound;
    double c = (a+b)/2.0;
    double d = (a+c)/2.0;
    double e = (c+b)/2.0;
    
    if (m_lowerBoundsInnerDimensions.size() > 0)
    {
        AdaptiveSimpson test;
        test.setFunction(m_integrand);
        test.setInterval(m_lowerBoundsInnerDimensions,m_upperBoundsInnerDimensions);
        test.setMaximumRecursionDepth(m_maxRecursionDepth);
        test.setPrecision(m_epsilon);
        m_evaluationPointsOuterDimensions.push_front(a);
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        first.fa = test.integrate();
        m_evaluationPointsOuterDimensions[0] = b;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        first.fb = test.integrate();
        m_evaluationPointsOuterDimensions[0] = c;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        first.fc = test.integrate();
        m_evaluationPointsOuterDimensions[0] = d;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        first.fd = test.integrate();
        m_evaluationPointsOuterDimensions[0] = e;
        test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
        first.fe = test.integrate();
    }
    else
    {
        m_evaluationPointsOuterDimensions.push_front(a);
        first.fa = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = b;
        first.fb = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = c;
        first.fc = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = d;
        first.fd = m_integrand(m_evaluationPointsOuterDimensions);
        m_evaluationPointsOuterDimensions[0] = e;
        first.fe = m_integrand(m_evaluationPointsOuterDimensions);
    }

    std::size_t size = first.fa.size();
    first.S.reserve(size);
    first.Sleft.reserve(size);
    first.Sright.reserve(size);
    double prefactor = (b-a)/12.0;
    for (std::size_t i = 0; i<size; i++)
    {
        first.S.emplace_back(2.0*prefactor*(first.fa[i]+4.0*first.fb[i]+first.fc[i]));
        first.Sleft.emplace_back(prefactor*(first.fa[i] + 4.0*first.fd[i] + first.fc[i]));
        first.Sright.emplace_back(prefactor*(first.fc[i] + 4.0*first.fe[i] + first.fb[i]));
        first.error = std::max(first.error,std::abs(15.0*(first.Sleft[i]+first.Sright[i]-first.S[i])));
    }
    first.lowerlimit=m_lowerBound;
    first.upperlimit=m_upperBound;
    first.epsilon = m_epsilon;
    
    std::priority_queue<helper,std::vector<helper>> myqueue;
    myqueue.emplace(first);
    
    helper element1,element2;    
    while (myqueue.size() < m_maxRecursionDepth )
    {
        auto mostError = myqueue.top();
        if (mostError.error < mostError.epsilon)
            break;
        
        //std::cout << mostError.error << " " << mostError.epsilon << std::endl;

        
        std::tie(element1,element2) = splitElement(mostError);
        myqueue.pop();
        myqueue.emplace(element1);
        myqueue.emplace(element2);
    }
    
    return sumPieces(myqueue);
}

AdaptiveSimpson::AdaptiveSimpson() : m_p{ new SimpsonImpl{} } {}

AdaptiveSimpson::AdaptiveSimpson(const AdaptiveSimpson& other) : m_p(other.m_p->clone()) {}

AdaptiveSimpson& AdaptiveSimpson::operator= (const AdaptiveSimpson& other)
{
    m_p = move(other.m_p->clone());
    return *this;
}

AdaptiveSimpson::AdaptiveSimpson(AdaptiveSimpson&& other)
{
    m_p = move(other.m_p);
    other.m_p = NULL;
}

AdaptiveSimpson& AdaptiveSimpson::operator= (AdaptiveSimpson&& other)
{
    if (m_p != other.m_p)
    {
        m_p = move(other.m_p);
        other.m_p = NULL;
    }
    return *this;
}

void AdaptiveSimpson::setFunction(const std::function< std::vector<double>(std::deque<double>& evaluationPoints)> & integrand)
{
    m_p->m_integrand = integrand;
}

void AdaptiveSimpson::setInterval(const std::vector<double>& lowerBounds, const std::vector<double>& upperBounds)
{
    assert(lowerBounds.size() == upperBounds.size());
    m_p->m_lowerBound = lowerBounds.back();
    m_p->m_upperBound = upperBounds.back();
    if (lowerBounds.size()>1)
    {
        m_p->m_lowerBoundsInnerDimensions = lowerBounds;
        m_p->m_lowerBoundsInnerDimensions.pop_back();
        
        m_p->m_upperBoundsInnerDimensions = upperBounds;
        m_p->m_upperBoundsInnerDimensions.pop_back();
    }
}

void AdaptiveSimpson::setPrecision(double epsilon)
{
    m_p->m_epsilon = epsilon;
}

void AdaptiveSimpson::setMaximumRecursionDepth(int maxRecursionDepth)
{
    m_p->m_maxRecursionDepth = maxRecursionDepth;
}
std::vector<double> AdaptiveSimpson::integrate()
{
    return m_p->integrate();
}

void AdaptiveSimpson::setAdditionalEvaluationPoints(const std::deque<double>& evaluationPoints)
{
    m_p->m_evaluationPointsOuterDimensions = evaluationPoints;
}

AdaptiveSimpson::~AdaptiveSimpson() {};

