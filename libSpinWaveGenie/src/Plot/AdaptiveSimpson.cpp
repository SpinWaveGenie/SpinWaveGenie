#include "AdaptiveSimpson.h"

#include <queue>
#include <algorithm>
#include <limits>

struct helper
{
  helper() : error(0.0){};
  std::vector<double> fa, fb, fc, fd, fe;
  std::vector<double> S, Sleft, Sright;
  double lowerlimit;
  double upperlimit;
  double epsilon;
  double error;
  bool operator<(const helper &rhs) const { return this->error / this->epsilon < rhs.error / rhs.epsilon; }
};

class AdaptiveSimpson::SimpsonImpl
{
public:
  SimpsonImpl() : m_epsilon(1.0e-5), m_maximumDivisions(1000){};
  std::vector<double> sumPieces(std::priority_queue<helper> &pieces);
  void createElement(helper &mostError, helper &element1);
  void splitElement(helper &mostError, helper &element1);
  std::vector<double> integrate();
  std::function<std::vector<double>(std::deque<double> &evaluationPoints)> m_integrand;
  double m_lowerBound, m_upperBound, m_epsilon;
  std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
  std::deque<double> m_evaluationPointsOuterDimensions;
  int m_maximumDivisions;
  std::unique_ptr<SimpsonImpl> clone();
};

std::unique_ptr<AdaptiveSimpson::SimpsonImpl> AdaptiveSimpson::SimpsonImpl::clone()
{
  return std::unique_ptr<AdaptiveSimpson::SimpsonImpl>(new SimpsonImpl(*this));
}

void AdaptiveSimpson::SimpsonImpl::createElement(helper &mostError,helper &element)
{
  double a1 = mostError.lowerlimit;
  double b1 = 0.5*(mostError.lowerlimit + mostError.upperlimit);
  double c1 = 0.5*(a1 + b1);
  double d1 = 0.5*(a1 + c1);
  double e1 = 0.5*(c1 + b1);
  double a2 = 0.5*(mostError.lowerlimit + mostError.upperlimit);
  double b2 = mostError.upperlimit;
  double c2 = 0.5*(a2 + b2);
  double d2 = 0.5*(a2 + c2);
  double e2 = 0.5*(c2 + b2);

  element.fa = std::move(mostError.fa);
  element.fb = mostError.fc;
  element.fc = std::move(mostError.fd);
  element.S = std::move(mostError.Sleft);
  
  mostError.fa = std::move(mostError.fc);
  mostError.fc = std::move(mostError.fe);
  mostError.S = std::move(mostError.Sright);

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    AdaptiveSimpson test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumDivisions(m_maximumDivisions);
    test.setPrecision(m_epsilon);
    m_evaluationPointsOuterDimensions[0] = d1;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fd = test.integrate();
    m_evaluationPointsOuterDimensions[0] = e1;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fe = test.integrate();
    m_evaluationPointsOuterDimensions[0] = d2;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    mostError.fd = test.integrate();
    m_evaluationPointsOuterDimensions[0] = e2;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    mostError.fe = test.integrate();
  }
  else
  {
    m_evaluationPointsOuterDimensions[0] = d1;
    element.fd = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = e1;
    element.fe = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = d2;
    mostError.fd = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = e2;
    mostError.fe = m_integrand(m_evaluationPointsOuterDimensions);
  }

  std::size_t size = element.fa.size();

  double prefactor = (b1 - a1) / 12.0;
  element.Sleft.reserve(size);
  element.Sright.reserve(size);
  
  mostError.Sleft.reserve(size);
  mostError.Sright.reserve(size);
  
  mostError.error = 0.0;
  for (std::size_t i = 0; i < size; i++)
  {
    element.Sleft.emplace_back(prefactor * (element.fa[i] + 4.0 * element.fd[i] + element.fc[i]));
    element.Sright.emplace_back(prefactor * (element.fc[i] + 4.0 * element.fe[i] + element.fb[i]));
    element.error = std::max(element.error, std::abs(15.0 * (element.Sleft[i] + element.Sright[i] - element.S[i])));
    
    mostError.Sleft.emplace_back(prefactor * (mostError.fa[i] + 4.0 * mostError.fd[i] + mostError.fc[i]));
    mostError.Sright.emplace_back(prefactor * (mostError.fc[i] + 4.0 * mostError.fe[i] + mostError.fb[i]));
    mostError.error = std::max(element.error, std::abs(15.0 * (mostError.Sleft[i] + mostError.Sright[i] - mostError.S[i])));
  }
  element.lowerlimit = a1;
  element.upperlimit = b1;
  element.epsilon = std::max(mostError.epsilon / sqrt(2.0), std::numeric_limits<double>::epsilon());
  
  mostError.lowerlimit = a2;
  mostError.upperlimit = b2;
  mostError.epsilon = std::max(mostError.epsilon / sqrt(2.0), std::numeric_limits<double>::epsilon());
}

void AdaptiveSimpson::SimpsonImpl::splitElement(helper &mostError,helper& element)
{
  this->createElement(mostError,element);
}

std::vector<double> AdaptiveSimpson::SimpsonImpl::sumPieces(std::priority_queue<helper> &pieces)
{
  std::size_t size = pieces.top().Sleft.size();
  std::vector<double> sum(size);
  while (pieces.size() > 0)
  {
    const helper& element = pieces.top();
    for (std::size_t i = 0; i < size; i++)
    {
      double S2 = element.Sleft[i] + element.Sright[i];
      sum[i] += S2 + (S2 - element.S[i]) / 15.0;
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
  double c = (a + b) / 2.0;
  double d = (a + c) / 2.0;
  double e = (c + b) / 2.0;

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    AdaptiveSimpson test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumDivisions(m_maximumDivisions);
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
  double prefactor = (b - a) / 12.0;
  for (std::size_t i = 0; i < size; i++)
  {
    first.S.emplace_back(2.0 * prefactor * (first.fa[i] + 4.0 * first.fb[i] + first.fc[i]));
    first.Sleft.emplace_back(prefactor * (first.fa[i] + 4.0 * first.fd[i] + first.fc[i]));
    first.Sright.emplace_back(prefactor * (first.fc[i] + 4.0 * first.fe[i] + first.fb[i]));
    first.error = std::max(first.error, std::abs(15.0 * (first.Sleft[i] + first.Sright[i] - first.S[i])));
  }
  first.lowerlimit = m_lowerBound;
  first.upperlimit = m_upperBound;
  first.epsilon = m_epsilon;

  std::priority_queue<helper, std::vector<helper>> myqueue;
  myqueue.emplace(first);

  while (myqueue.size() < m_maximumDivisions)
  {
    helper mostError(myqueue.top());
    if (mostError.error < mostError.epsilon)
      break;

    // std::cout << mostError.error << " " << mostError.epsilon << std::endl;
    helper element;
    splitElement(mostError,element);
    myqueue.pop();
    myqueue.emplace(std::move(mostError));
    myqueue.emplace(std::move(element));
  }

  return sumPieces(myqueue);
}

AdaptiveSimpson::AdaptiveSimpson() : m_p{new SimpsonImpl{}} {}

AdaptiveSimpson::AdaptiveSimpson(const AdaptiveSimpson &other) : m_p(other.m_p->clone()) {}

AdaptiveSimpson &AdaptiveSimpson::operator=(const AdaptiveSimpson &other)
{
  m_p = move(other.m_p->clone());
  return *this;
}

AdaptiveSimpson::AdaptiveSimpson(AdaptiveSimpson &&other)
{
  m_p = move(other.m_p);
  other.m_p = NULL;
}

AdaptiveSimpson &AdaptiveSimpson::operator=(AdaptiveSimpson &&other)
{
  if (m_p != other.m_p)
  {
    m_p = move(other.m_p);
    other.m_p = NULL;
  }
  return *this;
}

void
AdaptiveSimpson::setFunction(const std::function<std::vector<double>(std::deque<double> &evaluationPoints)> &integrand)
{
  m_p->m_integrand = integrand;
}

void AdaptiveSimpson::setInterval(const std::vector<double> &lowerBounds, const std::vector<double> &upperBounds)
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

void AdaptiveSimpson::setPrecision(double epsilon) { m_p->m_epsilon = epsilon; }

void AdaptiveSimpson::setMaximumDivisions(int maximumDivisions) { m_p->m_maximumDivisions = maximumDivisions; }
std::vector<double> AdaptiveSimpson::integrate() { return m_p->integrate(); }

void AdaptiveSimpson::setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints)
{
  m_p->m_evaluationPointsOuterDimensions = evaluationPoints;
}

AdaptiveSimpson::~AdaptiveSimpson(){};
