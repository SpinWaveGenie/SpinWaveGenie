#include "SpinWaveGenie/Plot/AdaptiveSimpson.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <queue>
#include <algorithm>
#include <limits>

struct helper
{
  helper() : lowerlimit(0.0), upperlimit(0.0), c(0.0), d(0.0), e(0.0), epsilon(1.0e-5), error(0.0){};
  double lowerlimit, upperlimit, c, d, e;
  std::vector<double> fa, fb, fc, fd, fe;
  std::vector<double> S, Sleft, Sright;
  double epsilon;
  double error;
  void resetBounds(double lowerlimit, double upperlimit);
  void initializeError();
  void updateError();
};

struct ComparePointers
{
  bool operator()(const std::shared_ptr<helper> &lhs, const std::shared_ptr<helper> &rhs) const
  {
    return lhs->error / lhs->epsilon < rhs->error / rhs->epsilon;
  }
};

class AdaptiveSimpson::SimpsonImpl
{
public:
  SimpsonImpl() : m_lowerBound(0.0), m_upperBound(0.0), m_epsilon(1.0e-5), m_maximumDivisions(1000){};
  std::vector<double> sumPieces(
      std::priority_queue<std::shared_ptr<helper>, std::vector<std::shared_ptr<helper>>, ComparePointers> &pieces);
  void createElement(const std::shared_ptr<helper> &mostError, const std::shared_ptr<helper> &element1);
  void splitElement(const std::shared_ptr<helper> &mostError, const std::shared_ptr<helper> &element1);
  std::vector<double> integrate();
  std::function<std::vector<double>(std::deque<double> &evaluationPoints)> m_integrand;
  double m_lowerBound, m_upperBound, m_epsilon;
  std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
  std::deque<double> m_evaluationPointsOuterDimensions;
  std::size_t m_maximumDivisions;
  std::unique_ptr<SimpsonImpl> clone();
};

std::unique_ptr<AdaptiveSimpson::SimpsonImpl> AdaptiveSimpson::SimpsonImpl::clone()
{
  return std::unique_ptr<AdaptiveSimpson::SimpsonImpl>(new SimpsonImpl(*this));
}

void helper::resetBounds(double lower, double upper)
{
  lowerlimit = lower;
  upperlimit = upper;
  c = 0.5 * (lowerlimit + upperlimit);
  d = 0.5 * (lowerlimit + c);
  e = 0.5 * (c + upperlimit);
}

void helper::updateError()
{
  error = 0.0;
  const std::size_t size = fa.size();
  double prefactor = (upperlimit - lowerlimit) / 12.0;
  Sleft.reserve(size);
  Sright.reserve(size);
  for (std::size_t i = 0; i < fa.size(); i++)
  {
    Sleft.emplace_back(prefactor * (fa[i] + 4.0 * fd[i] + fc[i]));
    Sright.emplace_back(prefactor * (fc[i] + 4.0 * fe[i] + fb[i]));
    error = std::max(error, std::abs(15.0 * (Sleft[i] + Sright[i] - S[i])));
  }
}

void helper::initializeError()
{
  const std::size_t size = fa.size();
  S.reserve(size);
  Sleft.reserve(size);
  Sright.reserve(size);
  double prefactor = (upperlimit - lowerlimit) / 12.0;
  for (std::size_t i = 0; i < size; i++)
  {
    double tmp = fa[i] + fc[i];
    S.emplace_back(2.0 * prefactor * (tmp + 4.0 * fb[i]));
    Sleft.emplace_back(prefactor * (tmp + 4.0 * fd[i]));
    Sright.emplace_back(prefactor * (fc[i] + 4.0 * fe[i] + fb[i]));
    error = std::max(error, std::abs(15.0 * (Sleft[i] + Sright[i] - S[i])));
  }
}

void AdaptiveSimpson::SimpsonImpl::createElement(const std::shared_ptr<helper> &mostError,
                                                 const std::shared_ptr<helper> &element)
{
  // element
  element->resetBounds(mostError->lowerlimit, 0.5 * (mostError->lowerlimit + mostError->upperlimit));
  element->epsilon = std::max(mostError->epsilon * M_SQRT1_2, std::numeric_limits<double>::epsilon());

  element->fa = std::move(mostError->fa);
  element->fb = mostError->fc;
  element->fc = std::move(mostError->fd);
  element->S = std::move(mostError->Sleft);

  // mostError
  mostError->resetBounds(0.5 * (mostError->lowerlimit + mostError->upperlimit), mostError->upperlimit);
  mostError->epsilon = std::max(mostError->epsilon * M_SQRT1_2, std::numeric_limits<double>::epsilon());

  mostError->fa = std::move(mostError->fc);
  mostError->fc = std::move(mostError->fe);
  mostError->S = std::move(mostError->Sright);

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    // both
    AdaptiveSimpson test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumDivisions(m_maximumDivisions);
    test.setPrecision(m_epsilon);
    // element
    m_evaluationPointsOuterDimensions[0] = element->d;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element->fd = test.integrate();
    m_evaluationPointsOuterDimensions[0] = element->e;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element->fe = test.integrate();
    // mostError
    m_evaluationPointsOuterDimensions[0] = mostError->d;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    mostError->fd = test.integrate();
    m_evaluationPointsOuterDimensions[0] = mostError->e;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    mostError->fe = test.integrate();
  }
  else
  {
    // element
    m_evaluationPointsOuterDimensions[0] = element->d;
    element->fd = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = element->e;
    element->fe = m_integrand(m_evaluationPointsOuterDimensions);
    // mostError
    m_evaluationPointsOuterDimensions[0] = mostError->d;
    mostError->fd = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = mostError->e;
    mostError->fe = m_integrand(m_evaluationPointsOuterDimensions);
  }

  element->updateError();
  mostError->updateError();
}

void AdaptiveSimpson::SimpsonImpl::splitElement(const std::shared_ptr<helper> &mostError,
                                                const std::shared_ptr<helper> &element)
{
  this->createElement(mostError, element);
}

std::vector<double> AdaptiveSimpson::SimpsonImpl::sumPieces(
    std::priority_queue<std::shared_ptr<helper>, std::vector<std::shared_ptr<helper>>, ComparePointers> &pieces)
{
  std::size_t size = pieces.top()->Sleft.size();
  std::vector<double> sum(size);
  double prefactor = 1.0 / 15.0;
  while (pieces.size() > 0)
  {
    const auto &element = pieces.top();
    for (std::size_t i = 0; i < size; i++)
    {
      double S2 = element->Sleft[i] + element->Sright[i];
      sum[i] += S2 + prefactor * (S2 - element->S[i]);
    }
    pieces.pop();
  }
  return sum;
}

std::vector<double> AdaptiveSimpson::SimpsonImpl::integrate()
{
  std::shared_ptr<helper> first = std::make_shared<helper>();
  first->resetBounds(m_lowerBound, m_upperBound);
  first->epsilon = m_epsilon;

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    AdaptiveSimpson test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumDivisions(m_maximumDivisions);
    test.setPrecision(m_epsilon);
    m_evaluationPointsOuterDimensions.push_front(first->lowerlimit);
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first->fa = test.integrate();
    m_evaluationPointsOuterDimensions[0] = first->upperlimit;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first->fb = test.integrate();
    m_evaluationPointsOuterDimensions[0] = first->c;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first->fc = test.integrate();
    m_evaluationPointsOuterDimensions[0] = first->d;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first->fd = test.integrate();
    m_evaluationPointsOuterDimensions[0] = first->e;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first->fe = test.integrate();
  }
  else
  {
    m_evaluationPointsOuterDimensions.push_front(first->lowerlimit);
    first->fa = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = first->upperlimit;
    first->fb = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = first->c;
    first->fc = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = first->d;
    first->fd = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = first->e;
    first->fe = m_integrand(m_evaluationPointsOuterDimensions);
  }

  first->initializeError();

  std::priority_queue<std::shared_ptr<helper>, std::vector<std::shared_ptr<helper>>, ComparePointers> myqueue(
      (ComparePointers()));
  myqueue.push(std::move(first));

  while (myqueue.size() < m_maximumDivisions)
  {
    std::shared_ptr<helper> mostError = myqueue.top();
    if (mostError->error < mostError->epsilon)
      break;
    std::shared_ptr<helper> element = std::make_shared<helper>();
    splitElement(mostError, element);
    myqueue.pop();
    myqueue.push(std::move(mostError));
    myqueue.push(std::move(element));
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

void AdaptiveSimpson::setMaximumDivisions(std::size_t maximumDivisions) { m_p->m_maximumDivisions = maximumDivisions; }
std::vector<double> AdaptiveSimpson::integrate() { return m_p->integrate(); }

void AdaptiveSimpson::setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints)
{
  m_p->m_evaluationPointsOuterDimensions = evaluationPoints;
}

AdaptiveSimpson::~AdaptiveSimpson(){};
