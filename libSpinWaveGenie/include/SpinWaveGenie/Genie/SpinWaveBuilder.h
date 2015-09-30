#ifndef __SpinWaveBuilder_H__
#define __SpinWaveBuilder_H__ 1

#include <string>
#include <vector>
#include "SpinWaveGenie/Memory.h"
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Interactions/Interaction.h"

namespace SpinWaveGenie
{

class SpinWaveBuilder
{
public:
  SpinWaveBuilder();
  SpinWaveBuilder(Cell &cellIn);
  void updateCell(Cell &cellIn);
  void addInteraction(std::unique_ptr<Interaction> in);
  void updateInteraction(std::string name, double value);
  double getEnergy();
  Eigen::VectorXcd getFirstOrderTerms();
  SpinWave createElement();

private:
  Cell cell;
  std::vector<std::unique_ptr<Interaction>> interactions;
};
}
#endif
