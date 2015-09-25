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
    class IntegrateNormalDistribution : public SpinWavePlot{
    public:
        IntegrateNormalDistribution()
        {
            m_Cell.setBasisVectors(1.0,1.0,1.0,90.0,90.0,90.0);
        };
        std::unique_ptr<SpinWavePlot> clone()
        {
          return std::unique_ptr<SpinWavePlot>(std::make_unique<IntegrateNormalDistribution>(*this));
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
        std::vector<double> getCut(double /*kx*/, double /*ky*/, double /*kz*/)
        {
            double frequency = 10.0;
            double FWHM = 1.0;
            std::vector<double> results;
            results.reserve(m_Energies.size());
            double ma = -4.0*log(2.0)/pow(FWHM,2);
            double factor = 2.0*sqrt(log(2.0))/(FWHM*sqrt(M_PI));
            for (auto it = m_Energies.begin(); it!=m_Energies.end(); ++it)
            {
                results.push_back(factor*exp(ma*pow(frequency-(*it),2)));

            }
            return results;
        };
        ~IntegrateNormalDistribution(){};
    private:
        SpinWaveGenie::Cell m_Cell;
        Energies m_Energies;
    };
}

BOOST_AUTO_TEST_CASE( NormalDistributionTest )
{
  std::unique_ptr<SpinWaveGenie::SpinWavePlot> res(std::make_unique<SpinWaveGenie::IntegrateNormalDistribution>());

    SpinWaveGenie::Energies centeredEnergies(10.0,10.0,2);
    std::unique_ptr<SpinWaveGenie::SpinWavePlot> cut(
        std::make_unique<SpinWaveGenie::IntegrateEnergy>(move(res), centeredEnergies, 3.0, 0.000001));
    std::vector<double> result = cut->getCut(0.0,0.0,1.0);
    BOOST_CHECK_CLOSE(result[0],1.0,1.0e-5);
    BOOST_CHECK_CLOSE(result[1],1.0,1.0e-5);
}



