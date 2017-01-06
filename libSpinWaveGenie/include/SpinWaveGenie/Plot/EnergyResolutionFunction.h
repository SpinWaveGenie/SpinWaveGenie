//
//  EnergyResolutionFunction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 2/5/14.
//
//

#ifndef __EnergyResolutionFunction__
#define __EnergyResolutionFunction__

#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Plot/OneDimensionalShapes.h"
#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include <algorithm>
#include <iostream>

namespace SpinWaveGenie
{

template <class T> class EnergyResolution : public SpinWavePlot
{
public:
  EnergyResolution() = default;
  EnergyResolution(const EnergyResolution &other);
  EnergyResolution &operator=(const EnergyResolution &other);
  EnergyResolution(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, const T &SWIn,
                   const Energies &energiesIn);
  std::vector<double> getCut(double kx, double ky, double kz) override;
  void setSpinWave(const T &SWIn);
  void setResolutionFunction(std::unique_ptr<OneDimensionalShapes> resolutionFunctionIn);
  const Cell &getCell() const override;
  void setEnergies(const Energies &energiesIn) override;
  const Energies &getEnergies() override;
  std::unique_ptr<SpinWavePlot> clone() override;

private:
  Energies energies;
  std::unique_ptr<OneDimensionalShapes> ResolutionFunction;
  T SW;
};

using EnergyResolutionFunction = EnergyResolution<SpinWaveGenie::SpinWave>;

template class EnergyResolution<SpinWave>;

template <class T>
EnergyResolution<T>::EnergyResolution(std::unique_ptr<OneDimensionalShapes> ResolutionFunctionIn, const T &SWIn,
                                      const Energies &energiesIn)
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
  {
    ResolutionFunction = std::move(other.ResolutionFunction->clone());
}
}

template <class T> EnergyResolution<T> &EnergyResolution<T>::operator=(const EnergyResolution &other)
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

template <class T> void EnergyResolution<T>::setSpinWave(const T &SWIn) { SW = SWIn; }

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

  for (const auto & point : points)
  {
    if (std::isfinite(point.frequency) && std::isfinite(point.intensity))
    {
      const double min = point.frequency + ResolutionFunction->getMinimumEnergy();
      const double max = point.frequency + ResolutionFunction->getMaximumEnergy();
      const std::size_t UpperBound = energies.getUpperBound(max);
      // std::cout << min << " " << energies.getLowerBound(min) << " " << max << " " << UpperBound << std::endl;
      for (std::size_t index = energies.getLowerBound(min); index != UpperBound; ++index)
      {
        fval[index] += point.intensity * ResolutionFunction->getFunction(point.frequency, energies[index]);
      }
    }
  }
  return fval;
}

template <class T> const Cell &EnergyResolution<T>::getCell() const { return SW.getCell(); }

template <class T> const Energies &EnergyResolution<T>::getEnergies() { return energies; }

template <class T> void EnergyResolution<T>::setEnergies(const Energies &energiesIn) { energies = energiesIn; }

template <class T> std::unique_ptr<SpinWavePlot> EnergyResolution<T>::clone()
{
  return memory::make_unique<EnergyResolution<T>>(*this);
}
}

#endif /* defined(__EnergyResolutionFunction__) */
