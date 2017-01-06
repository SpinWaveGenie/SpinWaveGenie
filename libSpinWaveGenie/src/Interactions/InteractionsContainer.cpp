#include "SpinWaveGenie/Interactions/InteractionsContainer.h"
#include "boost/format.hpp"
#include <algorithm> // std::sort
#include <cmath>

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

void InteractionsContainer::sort()
{
  std::sort(container.begin(), container.end(),
            [](const std::unique_ptr<Interaction> &a, const std::unique_ptr<Interaction> &b) { return *a < *b; });
}
std::ostream &operator<<(std::ostream &output, const SpinWaveGenie::InteractionsContainer &n)
{
  output << "  Name  sl1   sl2\n";
  for (const auto &result : n)
  {
    const auto sublattices = result.sublattices();
    output << result.getName() << " " << sublattices[0] << " " << sublattices[1] << "\n";
  }
  return output;
}
}
