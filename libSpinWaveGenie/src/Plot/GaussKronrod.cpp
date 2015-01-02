#include "GaussKronrod.h"
#include <queue>
#include <algorithm>
#include <limits>
#include <array>

class gaussianhelper
{
public:
  gaussianhelper() : error(0.0){};
  std::vector<double> S2;
  double lowerlimit;
  double upperlimit;
  double epsilon;
  double error;
  bool operator<(const gaussianhelper &rhs) const { return this->error / this->epsilon < rhs.error / rhs.epsilon; }
};

class GaussKronrod::GKImpl
{
public:
  GKImpl() : m_epsilon(1.0e-3), m_maxRecursionDepth(1000)
  {
    points =  {{-0.991455371120813,-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,-0.405845151377397,-0.207784955007898, 0.0,0.207784955007898,0.405845151377397,
      0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813}};
    gweights = {{0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870}};
    kweights = {{0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728,
      0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525, 0.104790010322250, 0.063092092629979, 0.022935322010529}};
  };
  std::array<double,7> gweights;
  std::array<double,15> points, kweights;
  std::vector<double> sumPieces(std::priority_queue<gaussianhelper> &pieces);
  gaussianhelper createElement(const gaussianhelper &mostError, bool first);
  std::pair<gaussianhelper, gaussianhelper> splitElement(const gaussianhelper &mostError);
  std::vector<double> integrate();
  std::function<std::vector<double>(std::deque<double> &evaluationPoints)> m_integrand;
  double m_lowerBound, m_upperBound, m_epsilon;
  std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
  std::deque<double> m_evaluationPointsOuterDimensions;
  int m_maxRecursionDepth;
  std::unique_ptr<GKImpl> clone();
};

std::unique_ptr<GaussKronrod::GKImpl> GaussKronrod::GKImpl::clone()
{
  return std::unique_ptr<GaussKronrod::GKImpl>(new GKImpl(*this));
}

struct getx
{
  double a,b;
  double operator()(double u)
  {
    return 0.5*(a+b - (b-a)*u);
  };
};

gaussianhelper GaussKronrod::GKImpl::createElement(const gaussianhelper &mostError, bool first)
{
  gaussianhelper element;
  getx pos;
  if (first)
  {
    pos.a = mostError.lowerlimit;
    pos.b = (mostError.lowerlimit + mostError.upperlimit) / 2.0;
  }
  else
  {
    pos.a = (mostError.lowerlimit + mostError.upperlimit) / 2.0;
    pos.b = mostError.upperlimit;
  }

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    GaussKronrod test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumRecursionDepth(m_maxRecursionDepth);
    test.setPrecision(m_epsilon);
  
    //i=0
    m_evaluationPointsOuterDimensions[0] = pos(points[0]);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    std::vector<double> ka = test.integrate();

    m_evaluationPointsOuterDimensions[0] = pos(points[1]);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    std::vector<double> ga = test.integrate();
    
    std::size_t size = ka.size();
    std::vector<double> S = std::vector<double>(size,0.0);
    element.S2 = std::vector<double>(size,0.0);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] += ga[j]*gweights[0];
      element.S2[j] += ka[j]*kweights[0] + ga[j]*kweights[1];
    }
    
    for (std::size_t i=1; i<7; i++)
    {
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i]);
      test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
      ka = test.integrate();
      
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i+1]);
      test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
      ga = test.integrate();
      
      for (std::size_t j = 0; j < size; j++)
      {
        S[j] += ga[j]*gweights[i];
        element.S2[j] += ka[j]*kweights[2*i] + ga[j]*kweights[2*i+1];
      }
    }
    
    //i = 7
    m_evaluationPointsOuterDimensions[0] = pos(points[14]);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    ka = test.integrate();
    double prefactor = 0.5*(pos.b-pos.a);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] *= prefactor;
      element.S2[j] += ka[j]*kweights[14];
      element.S2[j] *= prefactor;
      element.error = std::max(element.error,pow(200.0*std::abs(element.S2[j] - S[j]),1.5));
    }
  }
  else
  {
    
    //i=0
    m_evaluationPointsOuterDimensions[0] = pos(points[0]);
    std::vector<double> ka = m_integrand(m_evaluationPointsOuterDimensions);
    
    m_evaluationPointsOuterDimensions[0] = pos(points[1]);
    std::vector<double> ga = m_integrand(m_evaluationPointsOuterDimensions);
    
    std::size_t size = ka.size();
    std::vector<double> S = std::vector<double>(size,0.0);
    element.S2 = std::vector<double>(size,0.0);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] += ga[j]*gweights[0];
      element.S2[j] += ka[j]*kweights[0] + ga[j]*kweights[1];
    }
    
    for (std::size_t i=1; i<7; i++)
    {
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i]);
      ka = m_integrand(m_evaluationPointsOuterDimensions);
      
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i+1]);
      ga = m_integrand(m_evaluationPointsOuterDimensions);
      
      for (std::size_t j = 0; j < size; j++)
      {
        S[j] += ga[j]*gweights[i];
        element.S2[j] += ka[j]*kweights[2*i] + ga[j]*kweights[2*i+1];
      }
    }
    
    //i = 7
    m_evaluationPointsOuterDimensions[0] = pos(points[14]);
    ka = m_integrand(m_evaluationPointsOuterDimensions);
    double prefactor = 0.5*(pos.b-pos.a);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] *= prefactor;
      element.S2[j] += ka[j]*kweights[14];
      element.S2[j] *= prefactor;
      element.error = std::max(element.error,pow(200.0*std::abs(element.S2[j] - S[j]),1.5));
    }
  }
    
  element.lowerlimit = pos.a;
  element.upperlimit = pos.b;
  element.epsilon = std::max(mostError.epsilon * M_SQRT1_2, std::numeric_limits<double>::epsilon());
  return element;
}

std::pair<gaussianhelper, gaussianhelper> GaussKronrod::GKImpl::splitElement(const gaussianhelper &mostError)
{
  gaussianhelper element1 = this->createElement(mostError, true);
  gaussianhelper element2 = this->createElement(mostError, false);
  return std::make_pair(element1, element2);
}

std::vector<double> GaussKronrod::GKImpl::sumPieces(std::priority_queue<gaussianhelper> &pieces)
{
  std::size_t size = pieces.top().S2.size();
  std::vector<double> sum(size);
  while (pieces.size() > 0)
  {
    auto element = pieces.top();
    for (std::size_t i = 0; i < size; i++)
    {
      //double S2 = element.Sleft[i] + element.Sright[i];
      sum[i] += element.S2[i];
    }
    pieces.pop();
  }
  return sum;
}

std::vector<double> GaussKronrod::GKImpl::integrate()
{
  gaussianhelper element;
  
  getx pos;
  pos.a = m_lowerBound;
  pos.b = m_upperBound;
  
  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    GaussKronrod test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumRecursionDepth(m_maxRecursionDepth);
    test.setPrecision(m_epsilon);
    
    //i=0
    m_evaluationPointsOuterDimensions.push_front(pos(points[0]));
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    std::vector<double> ka = test.integrate();
    
    m_evaluationPointsOuterDimensions[0] = pos(points[1]);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    std::vector<double> ga = test.integrate();
    
    std::size_t size = ka.size();
    std::vector<double> S = std::vector<double>(size,0.0);
    element.S2 = std::vector<double>(size,0.0);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] += ga[j]*gweights[0];
      element.S2[j] += ka[j]*kweights[0] + ga[j]*kweights[1];
    }
    
    for (std::size_t i=1; i<7; i++)
    {
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i]);
      test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
      ka = test.integrate();
      
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i+1]);
      test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
      ga = test.integrate();
      
      for (std::size_t j = 0; j < size; j++)
      {
        S[j] += ga[j]*gweights[i];
        element.S2[j] += ka[j]*kweights[2*i] + ga[j]*kweights[2*i+1];
      }
    }
    
    //i = 7
    m_evaluationPointsOuterDimensions[0] = pos(points[14]);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    ka = test.integrate();
    double prefactor = 0.5*(pos.b-pos.a);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] *= prefactor;
      element.S2[j] += ka[j]*kweights[14];
      element.S2[j] *= prefactor;
      element.error = std::max(element.error,pow(200.0*std::abs(element.S2[j] - S[j]),1.5));
    }
  }
  else
  {
    
    //i=0
    m_evaluationPointsOuterDimensions.push_front(pos(points[0]));
    std::vector<double> ka = m_integrand(m_evaluationPointsOuterDimensions);
    
    m_evaluationPointsOuterDimensions[0] = pos(points[1]);
    std::vector<double> ga = m_integrand(m_evaluationPointsOuterDimensions);
    
    std::size_t size = ka.size();
    std::vector<double> S = std::vector<double>(size,0.0);
    element.S2 = std::vector<double>(size,0.0);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] += ga[j]*gweights[0];
      element.S2[j] += ka[j]*kweights[0] + ga[j]*kweights[1];
    }
    
    //i=1..6
    for (std::size_t i=1; i<7; i++)
    {
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i]);
      ka = m_integrand(m_evaluationPointsOuterDimensions);
      
      m_evaluationPointsOuterDimensions[0] = pos(points[2*i+1]);
      ga = m_integrand(m_evaluationPointsOuterDimensions);
      
      for (std::size_t j = 0; j < size; j++)
      {
        S[j] += ga[j]*gweights[i];
        element.S2[j] += ka[j]*kweights[2*i] + ga[j]*kweights[2*i+1];
      }
    }
    
    //i = 7
    m_evaluationPointsOuterDimensions[0] = pos(points[14]);
    ka = m_integrand(m_evaluationPointsOuterDimensions);
    double prefactor = 0.5*(pos.b-pos.a);
    for (std::size_t j = 0; j < size; j++)
    {
      S[j] *= prefactor;
      element.S2[j] += ka[j]*kweights[14];
      element.S2[j] *= prefactor;
      element.error = std::max(element.error,pow(200.0*std::abs(element.S2[j] - S[j]),1.5));
    }
  }
  element.lowerlimit = pos.a;
  element.upperlimit = pos.b;
  element.epsilon = m_epsilon;

  std::priority_queue<gaussianhelper, std::vector<gaussianhelper>> myqueue;
    myqueue.emplace(std::move(element));

  gaussianhelper element1, element2;
  while (myqueue.size() < m_maxRecursionDepth)
  {
    auto mostError = myqueue.top();
    if (mostError.error < mostError.epsilon)
      break;

    // std::cout << mostError.error << " " << mostError.epsilon << std::endl;

    std::tie(element1, element2) = splitElement(mostError);
    myqueue.pop();
    myqueue.emplace(std::move(element1));
    myqueue.emplace(std::move(element2));
  }

  return sumPieces(myqueue);
}

GaussKronrod::GaussKronrod() : m_p{new GKImpl{}} {}

GaussKronrod::GaussKronrod(const GaussKronrod &other) : m_p(other.m_p->clone()) {}

GaussKronrod &GaussKronrod::operator=(const GaussKronrod &other)
{
  m_p = move(other.m_p->clone());
  return *this;
}

GaussKronrod::GaussKronrod(GaussKronrod &&other)
{
  m_p = move(other.m_p);
  other.m_p = nullptr;
}

GaussKronrod &GaussKronrod::operator=(GaussKronrod &&other)
{
  if (m_p != other.m_p)
  {
    m_p = move(other.m_p);
    other.m_p = nullptr;
  }
  return *this;
}

void
GaussKronrod::setFunction(const std::function<std::vector<double>(std::deque<double> &evaluationPoints)> &integrand)
{
  m_p->m_integrand = integrand;
}

void GaussKronrod::setInterval(const std::vector<double> &lowerBounds, const std::vector<double> &upperBounds)
{
  assert(lowerBounds.size() == upperBounds.size());
  m_p->m_lowerBound = lowerBounds.back();
  m_p->m_upperBound = upperBounds.back();
  if (lowerBounds.size() > 1)
  {
    m_p->m_lowerBoundsInnerDimensions = lowerBounds;
    m_p->m_lowerBoundsInnerDimensions.pop_back();

    m_p->m_upperBoundsInnerDimensions = upperBounds;
    m_p->m_upperBoundsInnerDimensions.pop_back();
  }
}

void GaussKronrod::setPrecision(double epsilon) { m_p->m_epsilon = epsilon; }

void GaussKronrod::setMaximumRecursionDepth(int maxRecursionDepth) { m_p->m_maxRecursionDepth = maxRecursionDepth; }
std::vector<double> GaussKronrod::integrate() { return m_p->integrate(); }

void GaussKronrod::setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints)
{
  m_p->m_evaluationPointsOuterDimensions = evaluationPoints;
}

GaussKronrod::~GaussKronrod(){};
