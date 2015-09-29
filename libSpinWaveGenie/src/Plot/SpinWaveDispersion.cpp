#include <fstream>
#include "SpinWaveGenie/Plot/SpinWaveDispersion.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "External/ezRateProgressBar.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using namespace SpinWaveGenie;

SpinWaveDispersion::SpinWaveDispersion()
{
  PrintPosition = true;
  PrintFrequency = true;
  PrintIntensity = true;
}

void SpinWaveDispersion::setOptions(Options PrintOptions, bool Value)
{
  switch (PrintOptions)
  {
  case (Options::PrintPosition):
    PrintPosition = Value;
    break;
  case (Options::PrintFrequency):
    PrintFrequency = Value;
    break;
  case (Options::PrintIntensity):
    PrintIntensity = Value;
    break;
  }
}

void SpinWaveDispersion::setFilename(string name) { Filename = name; }

void SpinWaveDispersion::setGenie(const SpinWave &SW) { Genie = SW; }

void SpinWaveDispersion::setPoints(ThreeVectors<double> pts)
{
  for (auto it = pts.begin(); it != pts.end(); it++)
  {
    Kpoints.insert(it->get<0>(), it->get<1>(), it->get<2>());
  }
}

void SpinWaveDispersion::save()
{
  std::ofstream file(Filename);
  ez::ezRateProgressBar<std::size_t> p(Kpoints.size());
  p.units = "Q-points";
  p.start();
  for (auto it = Kpoints.begin(); it != Kpoints.end(); it++)
  {
    double x = it->get<0>();
    double y = it->get<1>();
    double z = it->get<2>();
    // cout << x << " " << y << " " << z << endl;

    if (PrintPosition)
    {
      file << x << " " << y << " " << z << " "; // << endl;
    }

    Genie.createMatrix(x, y, z);
    Genie.calculate();
    Results pts = Genie.getPoints();
    pts.sort();
    // pts.significantSolutions();

    if (PrintFrequency)
    {
      for (const auto & pt : pts)
      {
        file << pt.frequency << "  ";
      }
    }

    if (PrintIntensity)
    {
      for (const auto & pt : pts)
      {
        file << pt.intensity << " ";
      }
    }

    file << endl;
    p.update(std::distance(Kpoints.begin(), it));
  }
  p.update(Kpoints.size());
}
