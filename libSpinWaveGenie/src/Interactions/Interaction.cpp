#include "SpinWaveGenie/Interactions/Interaction.h"
using namespace std;

namespace SpinWaveGenie
{

bool Interaction::operator<(const Interaction &other) const
{
  auto sl1 = this->sublattices();
  auto sl2 = other.sublattices();
  for (std::size_t index = 0; index < sl1.size(); ++index)
  {
    if (sl1[index] < sl2[index])
      return true;
    else if (sl1[index] != sl2[index])
      return false;
  }
  return false;
}

}
