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
    BOOST_CHECK_CLOSE(result[0],6.0,1.0e-5);
    BOOST_CHECK_CLOSE(result[1],12.0,1.0e-5);
    BOOST_CHECK_CLOSE(result[2],18.0,1.0e-5);
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
                BOOST_CHECK_CLOSE(result[0],1.0,1.0e-5);
            }
            else
            {
                BOOST_CHECK_SMALL(result[0],1.0e-5);
            }
        }
    }
}
