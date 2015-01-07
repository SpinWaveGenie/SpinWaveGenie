#include "AdaptiveBoole.h"
#include <queue>
#include <algorithm>
#include <limits>

struct BooleHelper
{
  BooleHelper() : error(0.0){};
  std::vector<double> fa, fb, fc, fd, fe, ff, fg, fh, fi;
  std::vector<double> S, Sleft, Sright;
  double lowerlimit;
  double upperlimit;
  double epsilon;
  double error;
  bool operator<(const BooleHelper &rhs) const { return this->error / this->epsilon < rhs.error / rhs.epsilon; }
};

class AdaptiveBoole::BooleImpl
{
public:
  BooleImpl() : m_epsilon(1.0e-5), m_maximumDivisions(10){};
  std::vector<double> sumPieces(std::priority_queue<BooleHelper> &pieces);
  void createElement(const BooleHelper &mostError,BooleHelper &element, bool first);
  void splitElement(const BooleHelper &mostError,BooleHelper &element1,BooleHelper &element2);
  std::vector<double> integrate();
  std::function<std::vector<double>(std::deque<double> &evaluationPoints)> m_integrand;
  double m_lowerBound, m_upperBound, m_epsilon;
  std::vector<double> m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions;
  std::deque<double> m_evaluationPointsOuterDimensions;
  int m_maximumDivisions;
  std::unique_ptr<BooleImpl> clone();
};

std::unique_ptr<AdaptiveBoole::BooleImpl> AdaptiveBoole::BooleImpl::clone()
{
  return std::unique_ptr<AdaptiveBoole::BooleImpl>(new BooleImpl(*this));
}

void AdaptiveBoole::BooleImpl::createElement(const BooleHelper &mostError, BooleHelper& element, bool first)
{
  double a, b;
  if (first)
  {
    a = mostError.lowerlimit;
    b = 0.5*(mostError.lowerlimit + mostError.upperlimit);
  }
  else
  {
    a = 0.5*(mostError.lowerlimit + mostError.upperlimit);
    b = mostError.upperlimit;
  }
  double tmp = 0.125*(b - a);
  double c = a + tmp;
  double e = a + 3.0 * tmp;
  double g = a + 5.0 * tmp;
  double i = a + 7.0 * tmp;
  if (first)
  {
    element.fa = std::move(mostError.fa);
    element.fb = mostError.ff;
    element.fd = std::move(mostError.fc);
    element.ff = std::move(mostError.fd);
    element.fh = std::move(mostError.fe);
    element.S = std::move(mostError.Sleft);
  }
  else
  {
    element.fa = std::move(mostError.ff);
    element.fb = std::move(mostError.fb);
    element.fd = std::move(mostError.fg);
    element.ff = std::move(mostError.fh);
    element.fh = std::move(mostError.fi);
    element.S = std::move(mostError.Sright);
  }

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    AdaptiveBoole test;
    test.setFunction(m_integrand);
    test.setInterval(m_lowerBoundsInnerDimensions, m_upperBoundsInnerDimensions);
    test.setMaximumDivisions(m_maximumDivisions);
    test.setPrecision(m_epsilon);
    m_evaluationPointsOuterDimensions[0] = c;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fc = test.integrate();
    m_evaluationPointsOuterDimensions[0] = e;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fe = test.integrate();
    m_evaluationPointsOuterDimensions[0] = g;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fg = test.integrate();
    m_evaluationPointsOuterDimensions[0] = i;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    element.fi = test.integrate();
  }
  else
  {
    m_evaluationPointsOuterDimensions[0] = c;
    element.fc = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = e;
    element.fe = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = g;
    element.fg = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = i;
    element.fi = m_integrand(m_evaluationPointsOuterDimensions);
  }

  std::size_t size = element.fa.size();

  double prefactor = (b - a) / 180.0;
  element.Sleft.reserve(size);
  element.Sright.reserve(size);
  for (std::size_t i = 0; i < size; i++)
  {
    element.Sleft.emplace_back(prefactor * (7.0 * element.fa[i] + 32.0 * element.fc[i] + 12.0 * element.fd[i] +
                                            32.0 * element.fe[i] + 7.0 * element.ff[i]));
    element.Sright.emplace_back(prefactor * (7.0 * element.ff[i] + 32.0 * element.fg[i] + 12.0 * element.fh[i] +
                                             32.0 * element.fi[i] + 7.0 * element.fb[i]));
    element.error = std::max(element.error, std::abs(63.0 * (element.Sleft[i] + element.Sright[i] - element.S[i])));
  }
  element.lowerlimit = a;
  element.upperlimit = b;
  element.epsilon = std::max(mostError.epsilon * M_SQRT1_2, std::numeric_limits<double>::epsilon());
  // std::cout << element.epsilon << " " << element.error << std::endl;
}

void AdaptiveBoole::BooleImpl::splitElement(const BooleHelper &mostError, BooleHelper &element1, BooleHelper &element2)
{
  this->createElement(mostError, element1,true);
  this->createElement(mostError, element2,false);
  //return std::make_pair(element1, element2);
}

std::vector<double> AdaptiveBoole::BooleImpl::sumPieces(std::priority_queue<BooleHelper> &pieces)
{
  std::size_t size = pieces.top().Sleft.size();
  std::vector<double> sum(size);
  while (pieces.size() > 0)
  {
    const BooleHelper& element = pieces.top();
    for (std::size_t i = 0; i < size; i++)
    {
      double S2 = element.Sleft[i] + element.Sright[i];
      sum[i] += S2 + (S2 - element.S[i]) / 63.0;
    }
    pieces.pop();
  }
  return sum;
}

std::vector<double> AdaptiveBoole::BooleImpl::integrate()
{
  BooleHelper first;

  double a = m_lowerBound;
  double b = m_upperBound;

  double tmp = 0.125*(b - a);
  double c = a + tmp;
  double d = a + 2.0 * tmp;
  double e = a + 3.0 * tmp;
  double f = a + 4.0 * tmp;
  double g = a + 5.0 * tmp;
  double h = a + 6.0 * tmp;
  double i = a + 7.0 * tmp;
  //std::cout << b << " " << a + 8.0 * tmp << std::endl;
  

  if (m_lowerBoundsInnerDimensions.size() > 0)
  {
    AdaptiveBoole test;
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
    m_evaluationPointsOuterDimensions[0] = f;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first.ff = test.integrate();
    m_evaluationPointsOuterDimensions[0] = g;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first.fg = test.integrate();
    m_evaluationPointsOuterDimensions[0] = h;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first.fh = test.integrate();
    m_evaluationPointsOuterDimensions[0] = i;
    test.setAdditionalEvaluationPoints(m_evaluationPointsOuterDimensions);
    first.fi = test.integrate();
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
    m_evaluationPointsOuterDimensions[0] = f;
    first.ff = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = g;
    first.fg = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = h;
    first.fh = m_integrand(m_evaluationPointsOuterDimensions);
    m_evaluationPointsOuterDimensions[0] = i;
    first.fi = m_integrand(m_evaluationPointsOuterDimensions);
  }

  std::size_t size = first.fa.size();

  double prefactor = (b - a) / 180.0;
  first.S.reserve(size);
  first.Sleft.reserve(size);
  first.Sright.reserve(size);
  for (std::size_t i = 0; i < size; i++)
  {
    first.S.emplace_back(2.0 * prefactor * (7.0 * first.fa[i] + 32.0 * first.fd[i] + 12.0 * first.ff[i] +
                                            32.0 * first.fh[i] + 7.0 * first.fb[i]));
    first.Sleft.emplace_back(prefactor * (7.0 * first.fa[i] + 32.0 * first.fc[i] + 12.0 * first.fd[i] +
                                          32.0 * first.fe[i] + 7.0 * first.ff[i]));
    first.Sright.emplace_back(prefactor * (7.0 * first.ff[i] + 32.0 * first.fg[i] + 12.0 * first.fh[i] +
                                           32.0 * first.fi[i] + 7.0 * first.fb[i]));
    first.error = std::max(first.error, std::abs(63.0 * (first.Sleft[i] + first.Sright[i] - first.S[i])));
  }

  first.lowerlimit = m_lowerBound;
  first.upperlimit = m_upperBound;
  first.epsilon = m_epsilon;

  std::priority_queue<BooleHelper> myqueue;
  myqueue.emplace(first);

  while (myqueue.size() < m_maximumDivisions)
  {
    const BooleHelper& mostError = myqueue.top();
    if (mostError.error < mostError.epsilon)
      break;

    // std::cout << mostError.error << " " << mostError.epsilon << std::endl;

    // std::cout << "mostError: " << mostError.S[0] << " " << mostError.Sleft[0] << " " << mostError.Sright[0] << " "
    // << mostError.error << std::endl;
    BooleHelper element1, element2;
    splitElement(mostError,element1,element2);
    // std::cout << "element1: " << element1.S[0] << " " << element1.Sleft[0] << " " << element1.Sright[0] << " "  <<
    // element1.error << " "  << element1.epsilon << " " << myqueue.size() << std::endl;
    // std::cout << "element2: " << element2.S[0] << " " << element2.Sleft[0] << " " << element2.Sright[0] << " "  <<
    // element2.error << " " << element1.epsilon << std::endl;

    myqueue.pop();
    myqueue.emplace(std::move(element1));
    myqueue.emplace(std::move(element2));
  }

  return sumPieces(myqueue);
}

AdaptiveBoole::AdaptiveBoole() : m_p{new BooleImpl{}} {}

AdaptiveBoole::AdaptiveBoole(const AdaptiveBoole &other) : m_p(other.m_p->clone()) {}

AdaptiveBoole &AdaptiveBoole::operator=(const AdaptiveBoole &other)
{
  m_p = move(other.m_p->clone());
  return *this;
}

AdaptiveBoole::AdaptiveBoole(AdaptiveBoole &&other)
{
  m_p = move(other.m_p);
  other.m_p = NULL;
}

AdaptiveBoole &AdaptiveBoole::operator=(AdaptiveBoole &&other)
{
  if (m_p != other.m_p)
  {
    m_p = move(other.m_p);
    other.m_p = NULL;
  }
  return *this;
}

void
AdaptiveBoole::setFunction(const std::function<std::vector<double>(std::deque<double> &evaluationPoints)> &integrand)
{
  m_p->m_integrand = integrand;
}

void AdaptiveBoole::setInterval(const std::vector<double> &lowerBounds, const std::vector<double> &upperBounds)
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

void AdaptiveBoole::setPrecision(double epsilon) { m_p->m_epsilon = epsilon; }

void AdaptiveBoole::setMaximumDivisions(int maximumDivisions) { m_p->m_maximumDivisions = maximumDivisions; }
std::vector<double> AdaptiveBoole::integrate() { return m_p->integrate(); }

void AdaptiveBoole::setAdditionalEvaluationPoints(const std::deque<double> &evaluationPoints)
{
  m_p->m_evaluationPointsOuterDimensions = evaluationPoints;
}

AdaptiveBoole::~AdaptiveBoole(){};
