#define BOOST_TEST_MODULE FormFactorTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <vector>
#include <map>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Genie/MagneticFormFactor.h"

using namespace SpinWaveGenie;

class TestCoefficients : public MagneticFormFactor {
public:
    std::map<std::string,double> testCoefficients()
    {
        std::map<std::string,double> badCoefficients;
        for( auto it = coefficients.begin(); it!=coefficients.end(); it++)
        {
            std::vector<double> F = it->second;
            double test = F[0] + F[2] + F[4] + F[6];
            if(fabs(test-1.0)>5.0e-3)
            {
                badCoefficients.insert(std::pair<std::string,double>(it->first,test));
            }
        }
        return badCoefficients;
    }
};

BOOST_AUTO_TEST_CASE( CheckCoefficients )
{
    TestCoefficients test;
    std::map<std::string,double> badCoefficients;
    badCoefficients = test.testCoefficients();
    BOOST_CHECK_EQUAL(badCoefficients.size(), 0);
}

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    MagneticFormFactor FormFactor;
    BOOST_CHECK_THROW(FormFactor.getFormFactor(0.0,0.0,0.0),std::runtime_error);
}


BOOST_AUTO_TEST_CASE( AlternateConstructor )
{
    MagneticFormFactor FormFactor("NONE");
    
    for (int n = 0; n < 101; ++n)
    {
        BOOST_CHECK_CLOSE(FormFactor.getFormFactor(1.0,1.0,1.0),1.0,1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE( setType )
{
    MagneticFormFactor FormFactor;
    
    FormFactor.setType("FE3");
    
    for (int n = 0; n < 101; ++n)
    {
        double x = 0.0; double y = 50.0; double z = 0.0;
        
        std::vector<double> FE3 = { 0.397200, 13.244200,  0.629500,  4.903400, -0.031400,  0.349600,  0.004400};
        double s2 = (pow(x,2) + pow(y,2) + pow(z,2))/(16.0*M_PI*M_PI);
        double f_Q = FE3[0]*exp(-1.0*FE3[1]*s2) + FE3[2]*exp(-1.0*FE3[3]*s2) + FE3[4]*exp(-1.0*FE3[5]*s2) + FE3[6];
        
        BOOST_CHECK_CLOSE(FormFactor.getFormFactor(x,y,z),f_Q,1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE( setMultipleTypes )
{
    MagneticFormFactor FormFactor;
    
    std::vector<std::string> types = {"MN2","CO2"};
    std::vector<double> weights = {6.0,4.0};
        
    FormFactor.setType(types,weights);
        
    for (int n = 0; n < 101; ++n)
    {
        double x = 0.0; double y = 0.0; double z = n/50.0;
            
        std::vector<double> MN2 = { 0.422000, 17.684000,  0.594800,  6.005000,  0.004300, -0.609000, -0.021900};
        std::vector<double> CO2 = { 0.433200, 14.355300,  0.585700,  4.607700, -0.038200,  0.133800,  0.017900};

        double s2 = (pow(x,2) + pow(y,2) + pow(z,2))/(16.0*M_PI*M_PI);
        double MN_Q = MN2[0]*exp(-1.0*MN2[1]*s2) + MN2[2]*exp(-1.0*MN2[3]*s2) + MN2[4]*exp(-1.0*MN2[5]*s2) + MN2[6];
        double CO_Q = CO2[0]*exp(-1.0*CO2[1]*s2) + CO2[2]*exp(-1.0*CO2[3]*s2) + CO2[4]*exp(-1.0*CO2[5]*s2) + CO2[6];
       
        BOOST_CHECK_CLOSE(FormFactor.getFormFactor(x,y,z),0.6*MN_Q+0.4*CO_Q,1.0e-8);
    }
}

BOOST_AUTO_TEST_CASE( IncorrectNumberOfTypes )
{
    MagneticFormFactor FormFactor;
    
    std::vector<std::string> types = {"MN2","CO2"};
    std::vector<double> weights = {6.0,4.0,2.0};
    
    BOOST_CHECK_THROW(FormFactor.setType(types,weights),std::runtime_error);
}






