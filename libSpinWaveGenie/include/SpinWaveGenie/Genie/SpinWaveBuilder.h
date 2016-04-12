#ifndef __SpinWaveBuilder_H__
#define __SpinWaveBuilder_H__ 1

#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Interactions/InteractionsContainer.h"
#include "SpinWaveGenie/Memory.h"
#include <string>
#include <vector>

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
  InteractionsContainer interactions;
};
}
#endif
