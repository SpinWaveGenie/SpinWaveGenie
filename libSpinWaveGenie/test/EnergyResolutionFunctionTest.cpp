#define BOOST_TEST_MODULE EnergyResolutionFunctionTest
#define BOOST_TEST_MAIN
#include <iostream>
#include <numeric>
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Plot/Plot.h"

using namespace SpinWaveGenie;


class SimpleSpinWave
{
public:
  SimpleSpinWave() = default;
  SimpleSpinWave(Results input) { m_Results = input; };
  void createMatrix(double /*KX*/, double /*KY*/, double /*KZ*/){};
  void calculate(){};
  Results getPoints() { return m_Results; };
  const Cell &getCell() const { return m_Cell; };

private:
    Results m_Results;
    Cell m_Cell;
};

void runTest(std::unique_ptr<OneDimensionalShapes> resolutionFunction)
{
    Point pt;
    pt.frequency = 30.0;
    pt.intensity = 1.0;
    Results result;
    result.insert(pt);
    
    SimpleSpinWave SW(result);
    
    double max = resolutionFunction->getMaximumEnergy() + pt.frequency;
    double min = resolutionFunction->getMinimumEnergy() + pt.frequency;
    
    Energies energies(0.0, 60.0, 61);
    EnergyResolution<SimpleSpinWave> res(move(resolutionFunction),SW,energies);
    std::vector<double> testme = res.getCut(0.0, 0.0, 0.0);
    
    //check FWHM
    BOOST_CHECK_CLOSE(testme[25]/testme[30],0.5,0.01);
    
    //check lower bound;
    auto min_zeros = std::find_if(testme.begin(), testme.end(),std::bind(std::greater<>(),std::placeholders::_1,1.0e-3));
    BOOST_CHECK_EQUAL(std::distance(testme.begin(), min_zeros), static_cast<long>(energies.getLowerBound(min)));

    double shouldBeZero = std::accumulate(testme.begin(),min_zeros,0.0);
    BOOST_CHECK_CLOSE(shouldBeZero,0.0,1.0e-15);
    
    //check upper bound
    auto max_zeros = std::find_if_not(min_zeros+1, testme.end(),std::bind(std::greater<>(),std::placeholders::_1,1.0e-3));
    BOOST_CHECK_EQUAL(std::distance(testme.begin(), max_zeros), static_cast<long>(energies.getUpperBound(max)));

    shouldBeZero = std::accumulate(max_zeros,testme.end(),0.0);
    BOOST_CHECK_CLOSE(shouldBeZero,0.0,1.0e-15);
}


BOOST_AUTO_TEST_CASE( GaussianFunction )
{
    OneDimensionalFactory factory;
    auto gaussian = factory.getGaussian(10.0,5.0e-3);
    runTest(move(gaussian));
}

BOOST_AUTO_TEST_CASE( LorentzianFunction )
{
    OneDimensionalFactory factory;
    auto lorentzian = factory.getLorentzian(10.0,5.0e-3);
    runTest(move(lorentzian));
}

BOOST_AUTO_TEST_CASE(UpdatedGaussianFunction)
{
  OneDimensionalFactory factory;
  auto gaussian = factory.getGaussian(20.0, 1.0e-1);
  auto gaussian_ptr = dynamic_cast<OneDimensionalGaussian *>(gaussian.get());
  BOOST_CHECK(gaussian_ptr != nullptr);
  // fix coverity issue
  if (gaussian_ptr)
  {
    gaussian_ptr->setFWHM(10.0);
    gaussian_ptr->setTolerance(5.0e-3);
    runTest(move(gaussian));
  }
}

BOOST_AUTO_TEST_CASE(UpdatedLorentzianFunction)
{
  OneDimensionalFactory factory;
  auto lorentzian = factory.getLorentzian(20.0, 1.0e-1);
  auto lorentzian_ptr = dynamic_cast<OneDimensionalLorentzian *>(lorentzian.get());
  BOOST_CHECK(lorentzian_ptr != nullptr);
  // fix coverity issue
  if (lorentzian_ptr)
  {
    lorentzian_ptr->setFWHM(10.0);
    lorentzian_ptr->setTolerance(5.0e-3);
    runTest(move(lorentzian));
  }
}
