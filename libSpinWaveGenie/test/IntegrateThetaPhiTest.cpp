#define BOOST_TEST_MODULE IntegrateThetaPhiTest
#define BOOST_TEST_MAIN
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Plot/Plot.h"

using namespace SpinWaveGenie;

namespace SpinWaveGenie
{
class SphericalHarmonics : public SpinWavePlot{
public:
    SphericalHarmonics()
    {
        m_Cell.setBasisVectors(1.0,1.0,1.0,90.0,90.0,90.0);
        m_Energies = Energies(0.0,0.0,1);
    }
    std::unique_ptr<SpinWavePlot> clone()
    {
        return std::unique_ptr<SpinWavePlot>(new SphericalHarmonics(*this));
    };
    const Cell& getCell() const
    {
        return m_Cell;
    };
    const Energies& getEnergies()
    {
        return m_Energies;
    };
    void setEnergies(Energies energies)
    {
        m_Energies = energies;
    };
    std::vector<double> getCut(double kx, double ky, double kz)
    {
        double r = sqrt(kx*kx+ky*ky+kz*kz);
        double theta = acos(kz/r);
        double phi = atan2(ky,kx);
        double tmp1 = boost::math::spherical_harmonic_r(2, 0, theta,phi);
        double tmp2 = boost::math::spherical_harmonic_r(4, 0, theta,phi);
        return std::vector<double>(1,tmp1*tmp2);
    };
    ~SphericalHarmonics(){};
private:
    SpinWaveGenie::Cell m_Cell;
    Energies m_Energies;
};
}

BOOST_AUTO_TEST_CASE( SphericalHarmonicsTest )
{
    std::unique_ptr<SphericalHarmonics> res(new SphericalHarmonics());
    std::unique_ptr<SpinWavePlot> cut(new IntegrateThetaPhi(move(res),0.001));
    std::vector<double> result = cut->getCut(0.0,0.0,1.0);
    std::cout << result[0]*4.0*M_PI << std::endl; //result from IntegrateThetaPhi is divided by 4*M_PI
}

