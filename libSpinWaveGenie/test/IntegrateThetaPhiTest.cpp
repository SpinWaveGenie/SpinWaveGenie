#define BOOST_TEST_MODULE IntegrateThetaPhiTest
#define BOOST_TEST_MAIN
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Plot/Plot.h"

namespace SpinWaveGenie
{
    class ConstantFunction : public SpinWavePlot{
    public:
        ConstantFunction()
        {
            m_Cell.setBasisVectors(1.0,1.0,1.0,90.0,90.0,90.0);
            m_Energies = Energies(0.0,0.0,1);
        };
        std::unique_ptr<SpinWavePlot> clone() override
        {
            return memory::make_unique<ConstantFunction>(*this);
        };
        const Cell& getCell() const override
        {
            return m_Cell;
        };
        const Energies& getEnergies() override
        {
            return m_Energies;
        };
        void setEnergies(Energies energies) override
        {
            m_Energies = energies;
        };
        std::vector<double> getCut(double /*kx*/, double /*ky*/, double /*kz*/) override { return std::vector<double>(1, 1.0); };
        ~ConstantFunction(){};
    private:
        SpinWaveGenie::Cell m_Cell;
        Energies m_Energies;
    };
}

BOOST_AUTO_TEST_CASE( ConstantFunctionTest )
{
  std::unique_ptr<SpinWaveGenie::SpinWavePlot> res(
      SpinWaveGenie::memory::make_unique<SpinWaveGenie::ConstantFunction>());
  std::unique_ptr<SpinWaveGenie::SpinWavePlot> cut(
      SpinWaveGenie::memory::make_unique<SpinWaveGenie::IntegrateThetaPhi>(move(res), 1.0e-10));
    std::vector<double> result = cut->getCut(0.0,0.0,1.0);
    //result from IntegrateThetaPhi is divided by 4*M_PI
    BOOST_CHECK_CLOSE(result[0],1.0,1.0e-5);
    //check another radius.
    result = cut->getCut(0.0,0.0,2.0);
    BOOST_CHECK_CLOSE(result[0],1.0,1.0e-5);
}

namespace SpinWaveGenie
{
class SphericalHarmonics : public SpinWavePlot{
public:
    SphericalHarmonics(unsigned n1, unsigned n2)
    {
        m_n1 = n1;
        m_n2 = n2;
        m_Cell.setBasisVectors(1.0,1.0,1.0,90.0,90.0,90.0);
        m_Energies = Energies(0.0,0.0,1);
    };
    std::unique_ptr<SpinWavePlot> clone() override
    {
        return memory::make_unique<SphericalHarmonics>(*this);
    };
    const Cell& getCell() const override
    {
        return m_Cell;
    };
    const Energies& getEnergies() override
    {
        return m_Energies;
    };
    void setEnergies(Energies energies) override
    {
        m_Energies = energies;
    };
    std::vector<double> getCut(double kx, double ky, double kz) override
    {
        double r = sqrt(kx*kx+ky*ky+kz*kz);
        double theta = acos(kz/r);
        double phi = atan2(ky,kx);
        double tmp1 = boost::math::spherical_harmonic_r(m_n1, 0, theta,phi);
        double tmp2 = boost::math::spherical_harmonic_r(m_n2, 0, theta,phi);
        return std::vector<double>(1,tmp1*tmp2);
    };
private:
    SpinWaveGenie::Cell m_Cell;
    Energies m_Energies;
    unsigned m_n1,m_n2;
};
}

BOOST_AUTO_TEST_CASE( SphericalHarmonicsTest )
{
    for(unsigned n1=0;n1<4;++n1)
    {
        for(unsigned n2=0;n2<4;++n2)
        {
          std::unique_ptr<SpinWaveGenie::SphericalHarmonics> res(
              SpinWaveGenie::memory::make_unique<SpinWaveGenie::SphericalHarmonics>(n1, n2));
          std::unique_ptr<SpinWaveGenie::SpinWavePlot> cut(
              SpinWaveGenie::memory::make_unique<SpinWaveGenie::IntegrateThetaPhi>(move(res), 1.0e-12));
            std::vector<double> result = cut->getCut(0.0,0.0,1.0);
            if (n1 == n2)
            {
                //result from IntegrateThetaPhi is divided by 4*M_PI
                BOOST_CHECK_CLOSE(result[0]*4.0*M_PI,1.0,5.0e-3);
            }
            else
            {
                BOOST_CHECK_SMALL(result[0],1.0e-5);
            }
        }
    }
}



