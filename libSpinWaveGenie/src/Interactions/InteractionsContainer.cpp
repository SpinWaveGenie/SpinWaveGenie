#include "SpinWaveGenie/Interactions/InteractionsContainer.h"
#include "boost/format.hpp"
#include <algorithm> // std::sort
#include <cmath>

using std::vector;

namespace SpinWaveGenie
{
InteractionsContainer::InteractionsContainer(const InteractionsContainer &other)
{
  this->container.reserve(other.size());
  for (const auto &interaction : other.container)
  {
    this->container.push_back(interaction->clone());
  }
}

InteractionsContainer &InteractionsContainer::operator=(const InteractionsContainer &other)
{
  this->container.reserve(other.size());
  for (const auto &interaction : other.container)
  {
    this->container.push_back(interaction->clone());
  }
  return *this;
}

void InteractionsContainer::clear() { container.clear(); }

void InteractionsContainer::insert(std::unique_ptr<Interaction> value) { container.push_back(std::move(value)); }

std::ostream &operator<<(std::ostream &output, const SpinWaveGenie::InteractionsContainer &n)
{
  /*output << "  frequency  intensity\n";
  for (const auto &result : n)
  {
      output << boost::format("%9.5f %10.5f\n") % result.frequency % result.intensity;
  }*/
  return output;
}
}
