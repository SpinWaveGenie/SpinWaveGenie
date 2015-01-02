//
//  EnergyResolutionFunction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __EnergyResolutionFunction__
#define __EnergyResolutionFunction__

#include <iostream>
#include <memory>
#include <algorithm>
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Containers/Results.h"

namespace SpinWaveGenie
{

template <class T> class EnergyResolution : public SpinWavePlot
{
public:
  EnergyResolution(){};
  EnergyResolution(const EnergyResolution &other);
  EnergyResolution &operator=(EnergyResolution &other);
  EnergyResolution(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, T SWIn, Energies energies);
  std::vector<double> getCut(double kxIn, double kyIn, double kzIn);
  void setSpinWave(T SWIn);
  void setResolutionFunction(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn);
  const Cell &getCell() const;
  void setEnergies(Energies energies);
  const Energies &getEnergies();
  std::unique_ptr<SpinWavePlot> clone();
  ~EnergyResolution(){};

private:
  void calculateEnergies();
  std::size_t getBin(double Energy);
  Energies energies;
  std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
  T SW;
};

typedef EnergyResolution<SpinWave> EnergyResolutionFunction;

template class EnergyResolution<SpinWave>;

template <class T>
EnergyResolution<T>::EnergyResolution(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, T SWIn,
                                      Energies energiesIn)
{
  // std::cout << "Creating Energy Resolution Function" << std::endl;
  this->energies = energiesIn;
  // cout << "Energy Points " << EnergyPoints << endl;
  ResolutionFunction = std::move(ResolutionFunctionIn);
  SW = SWIn;
}

template <class T> EnergyResolution<T>::EnergyResolution(const EnergyResolution<T> &other)
{
  // std::cout << "Copying Energy Resolution Function" << std::endl;
  energies = other.energies;
  // cout << "Energy Points??? " << other.EnergyPoints << endl;
  // cout << "Energy Points??? " << EnergyPoints << endl;
  SW = other.SW;
  if (other.ResolutionFunction)
    ResolutionFunction = std::move(other.ResolutionFunction->clone());
}

template <class T> EnergyResolution<T> &EnergyResolution<T>::operator=(EnergyResolution &other)
{
  std::cout << "Copying Energy Resolution Function" << std::endl;
  energies = other.energies;
  SW = other.SW;
  ResolutionFunction = move(other.ResolutionFunction->clone());
  return *this;
}

template <class T>
void EnergyResolution<T>::setResolutionFunction(std::unique_ptr<OneDimensionalShapes> resolutionFunctionIn)
{
  ResolutionFunction = std::move(resolutionFunctionIn);
}

template <class T> void EnergyResolution<T>::setSpinWave(T SWIn) { SW = SWIn; }

template <class T> std::vector<double> EnergyResolution<T>::getCut(double kx, double ky, double kz)
{
  // cout << "Energy Points: " << EnergyPoints << endl;
  // cout << MinimumEnergy << " " << MaximumEnergy << endl;
  size_t EnergyPoints = energies.size();
  std::vector<double> fval(EnergyPoints, 0.0);

  SW.createMatrix(kx, ky, kz);
  SW.calculate();
  Results points = SW.getPoints();

  // points.significantSolutions();

  for (auto pt = points.begin(); pt != points.end(); pt++)
  {
    if (std::isnan(pt->frequency) || std::isnan(pt->intensity))
    {
      // cout << "found NaN: " << pt->frequency << " " << pt->intensity << endl;
    }
    else
    {
      double min = pt->frequency + ResolutionFunction->getMinimumEnergy();
      double max = pt->frequency + ResolutionFunction->getMaximumEnergy();
      size_t UpperBound = energies.getUpperBound(max);
      // std::cout << min << " " << energies.getLowerBound(min) << " " << max << " " << UpperBound << std::endl;
      for (size_t index = energies.getLowerBound(min); index != UpperBound; index++)
      {
        fval[index] += pt->intensity * ResolutionFunction->getFunction(pt->frequency, energies[index]);
      }
    }
  }
  return fval;
}

template <class T> const Cell &EnergyResolution<T>::getCell() const { return SW.getCell(); }

template <class T> const Energies &EnergyResolution<T>::getEnergies() { return energies; }

template <class T> void EnergyResolution<T>::setEnergies(Energies energiesIn) { energies = energiesIn; }

template <class T> std::unique_ptr<SpinWavePlot> EnergyResolution<T>::clone()
{
  return std::unique_ptr<SpinWavePlot>(new EnergyResolution<T>(*this));
}
}

#endif /* defined(__EnergyResolutionFunction__) */
