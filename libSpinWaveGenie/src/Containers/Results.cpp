//
//  Results.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/1/14.
//
//

#include "SpinWaveGenie/Containers/Results.h"
#include "boost/format.hpp"
#include <algorithm> // std::sort
#include <cmath>

using std::vector;

namespace SpinWaveGenie
{

bool evalues_equal(const Point &a, const Point &b)
{
  // remove eigenvalues that are equal
  double EPS = 1.0e-5;
  return std::abs(a.frequency - b.frequency) < EPS;
}

void Results::clear() { results.clear(); }

void Results::insert(Point value) { results.push_back(value); }

void Results::sort() { std::sort(results.begin(), results.end()); }

void Results::uniqueSolutions()
{
  double EPS = 1.0e-5;
  vector<Point> VI_unique;
  // Find unique eigenvalues
  VI_unique = results;
  std::sort(VI_unique.begin(), VI_unique.end());
  VI_unique.erase(unique(VI_unique.begin(), VI_unique.end(), evalues_equal), VI_unique.end());

  std::size_t NU = VI_unique.size();
  VI_unique.resize(NU);

  for (std::size_t i = 0; i < NU; i++)
  {
    VI_unique[i].intensity = 0.0;
  }

  for (const auto &elem : results)
  {
    std::size_t VP_pos = NU; // set position to a nonsense value
    for (std::size_t j = 0; j < NU; j++)
    {
      if (std::abs(elem.frequency - VI_unique[j].frequency) < EPS)
      {
        VP_pos = j;
        VI_unique[j].intensity += elem.intensity;
        break;
      }
    }
    if (VP_pos == NU)
    {
      std::cout << "error finding unique value" << std::endl;
    }
  }

  results = VI_unique;
}

void Results::significantSolutions(double ETS)
{
  vector<Point> VI_signif;
  for (size_t k = 0; k != results.size(); k++)
  {
    if (results[k].intensity > ETS)
    {
      VI_signif.push_back(results[k]);
    }
  }
  results = VI_signif;
}

std::ostream &operator<<(std::ostream &output, const SpinWaveGenie::Results &n)
{
  output << "  frequency  intensity\n";
  for (const auto &result : n)
  {
    output << boost::format("%9.5f %10.5f\n") % result.frequency % result.intensity;
  }
  return output;
}
} // namespace SpinWaveGenie
