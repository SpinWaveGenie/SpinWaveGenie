#define BOOST_TEST_MODULE AdaptiveSimpsonTest
#define BOOST_TEST_MAIN
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "SpinWaveGenie/Plot/AdaptiveSimpson.h"

std::vector<double> constantFunction(std::deque<double>& x)
{
    std::vector<double> result;
    result.reserve(x.size());
    for (std::size_t i = 0; i<x.size();i++)
    {
        result.push_back(i+1);
    }
    return result;
};

BOOST_AUTO_TEST_CASE( ConstantFunctionTest )
{
    std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(constantFunction,std::placeholders::_1);
    AdaptiveSimpson test;
    test.setFunction(f);
    std::vector<double> lower_bounds = {0.0,0.0,0.0};
    std::vector<double> upper_bounds = {1.0,2.0,3.0};
    test.setInterval(lower_bounds,upper_bounds);
    auto result = test.integrate();
    BOOST_CHECK_CLOSE(result[0],6.0,1.0e-3);
    BOOST_CHECK_CLOSE(result[1],12.0,1.0e-3);
    BOOST_CHECK_CLOSE(result[2],18.0,1.0e-3);
}

class SphericalHarmonics
{
public:
    SphericalHarmonics(unsigned n1, unsigned n2)
    {
        m_n1 = n1;
        m_n2 = n2;
    };
    std::vector<double> getPoint(std::deque<double>& angles)
    {
        double theta = angles[0];
        double phi = angles[1];
        double tmp1 = boost::math::spherical_harmonic_r(m_n1, 0, theta,phi);
        double tmp2 = boost::math::spherical_harmonic_r(m_n2, 0, theta,phi);
        return std::vector<double>(1,tmp1*tmp2*sin(theta));
    };
private:
    unsigned m_n1,m_n2;
};

BOOST_AUTO_TEST_CASE( SphericalHarmonicsTest )
{
    for(unsigned n1=0;n1<4;++n1)
    {
        for(unsigned n2=0;n2<4;++n2)
        {    
            SphericalHarmonics myClass(n1,n2);
            std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(&SphericalHarmonics::getPoint,&myClass,std::placeholders::_1);
            AdaptiveSimpson test;
            test.setFunction(f);
            std::vector<double> lower_bounds = {0.0,0.0};
            std::vector<double> upper_bounds = {M_PI,2.0*M_PI};
            test.setInterval(lower_bounds,upper_bounds);
            auto result = test.integrate();
    
	    if (n1 == n2)
	    {
                BOOST_CHECK_CLOSE(result[0],1.0,1.0e-3);
            }
            else
            {
                BOOST_CHECK_SMALL(result[0],1.0e-3);
            }
        }
    }
}

/* Simple product function */
//double f0 (unsigned dim, const double *x, void *params)
std::vector<double> productFunction(std::deque<double>& x)
{
  std::vector<double> prod(1,1.0);
  unsigned int i;
  for (i = 0; i < x.size(); ++i)
    prod[0] *= 2.0 * x[i];
  return prod;
}

/* Gaussian centered at 1/2. */
std::vector<double> gaussianFunction(std::deque<double>& x)
{
  std::size_t dim = x.size();
  double a = 0.1;
  double sum = 0.0;
  for ( unsigned int i = 0; i < dim; i++)
  {
    double dx = x[i] - 0.5;
    sum += dx * dx;
  }
  std::vector<double> result(1);
  result[0] = std::pow (M_2_SQRTPI / (2.0 * a),dim) * exp (-sum / (a * a));
  return result;
}

BOOST_AUTO_TEST_CASE( GaussianFunctionTest )
{
  std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(gaussianFunction,std::placeholders::_1);
  AdaptiveSimpson test;
  test.setFunction(f);
  std::vector<double> lower_bounds = {0.0,0.0};
  std::vector<double> upper_bounds = {1.0,1.0};
  test.setInterval(lower_bounds,upper_bounds);
  auto result = test.integrate();
  BOOST_CHECK_CLOSE(result[0],1.0,1.0e-3);
}

/* double gaussian */
std::vector<double> doubleGaussianFunction(std::deque<double>& x)
{
  std::size_t dim = x.size();
  double a = 0.1;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for ( unsigned int i = 0; i < dim; i++)
  {
    double dx1 = x[i] - 1.0 / 3.0;
    double dx2 = x[i] - 2.0 / 3.0;
    sum1 += dx1 * dx1;
    sum2 += dx2 * dx2;
  }
  std::vector<double> result(1);
  result[0] = 0.5 * std::pow (M_2_SQRTPI / (2.0 * a), dim) * (std::exp (-sum1 / (a * a)) + std::exp (-sum2 / (a * a)));
  return result;
}

//also tests move constructor
BOOST_AUTO_TEST_CASE( DoubleGaussianFunctionTest )
{
  std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(doubleGaussianFunction,std::placeholders::_1);
  AdaptiveSimpson test;
  test.setFunction(f);
  std::vector<double> lower_bounds = {0.0,0.0};
  std::vector<double> upper_bounds = {1.0,1.0};
  test.setInterval(lower_bounds,upper_bounds);
  auto result1 = test.integrate();
  BOOST_CHECK_CLOSE(result1[0],1.0,1.0e-3);
  AdaptiveSimpson moveTest(std::move(test));
  auto result2 = moveTest.integrate();
  BOOST_CHECK_CLOSE(result2[0],1.0,1.0e-3);
}

/* Tsuda's example */
std::vector<double> tsudasFunction(std::deque<double>& x)
{
  std::size_t dim = x.size();
  double c = (1.0+sqrt (10.0))/9.0;
  std::vector<double> result(1,1.0);
  unsigned int i;
  for (i = 0; i < dim; i++)
    result[0] *= c / (c + 1) * pow((c + 1) / (c + x[i]), 2.0);
  return result;
}

//also tests copy constructor
BOOST_AUTO_TEST_CASE( TsudasFunctionTest )
{
  std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(tsudasFunction,std::placeholders::_1);
  AdaptiveSimpson test;
  test.setFunction(f);
  std::vector<double> lower_bounds = {0.0,0.0};
  std::vector<double> upper_bounds = {1.0,1.0};
  test.setInterval(lower_bounds,upper_bounds);
  auto result1 = test.integrate();
  BOOST_CHECK_CLOSE(result1[0],1.0,1.0e-3);
  AdaptiveSimpson copyTest(test);
  auto result2 = copyTest.integrate();
  BOOST_CHECK_CLOSE(result2[0],1.0,1.0e-3);
}

std::vector<double> complexFunction(std::deque<double>& x)
{
  std::size_t dim = x.size();
  double p = 1.0 / dim;
  double prod = pow(1 + p, dim);
  unsigned int i;
  for (i = 0; i < dim; i++)
    prod *= pow(x[i], p);
  return std::vector<double>(1,prod);
};

//also tests copy assignment operator
BOOST_AUTO_TEST_CASE( ComplexFunctionTest )
{
  std::function< std::vector<double>(std::deque<double>& x)> f = std::bind<std::vector<double> >(complexFunction,std::placeholders::_1);
  AdaptiveSimpson test;
  test.setFunction(f);
  std::vector<double> lower_bounds = {0.0,0.0,0.0};
  std::vector<double> upper_bounds = {1.0,1.0,1.0};
  test.setInterval(lower_bounds,upper_bounds);
  AdaptiveSimpson copyTest = test;
  auto result1 = test.integrate();
  BOOST_CHECK_CLOSE(result1[0],1.0,1.0e-3);
  auto result2 = copyTest.integrate();
  BOOST_CHECK_CLOSE(result2[0],1.0,1.0e-3);
}
