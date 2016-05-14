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

void SpinWaveDispersion::setFilename(const std::string &name) { Filename = name; }

void SpinWaveDispersion::setGenie(const SpinWave &SW) { Genie = SW; }

void SpinWaveDispersion::setPoints(const ThreeVectors<double> &points)
{
  for (const auto &pt : points)
  {
    Kpoints.insert(pt[0], pt[1], pt[2]);
  }
}

void SpinWaveDispersion::save()
{
  std::ofstream file(Filename);
  ez::ezRateProgressBar<std::size_t> p(Kpoints.size());
  p.units = "Q-points";
  p.start();
  for (const auto &kpt : Kpoints)
  {
    double x = kpt[0];
    double y = kpt[1];
    double z = kpt[2];
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

    file << "\n";
    ++p;
  }
  p.update(Kpoints.size());
}
