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
  //! Default Constructor
  SpinWaveBuilder();
  //! Construct a SpinWaveBuilder object
  //! \param cellIn Unit cell for a given simulation.
  SpinWaveBuilder(const Cell &cellIn);
  //! Update the unit cell inside the SpinWaveBuilder object.
  //! \param cellIn New unit cell.
  void updateCell(const Cell &cellIn);
  //! Add a magnetic interaction to the model using move semantics (avoids a copy)
  //! \param in magnetic interaction being added to the model builder.
  void addInteraction(std::unique_ptr<Interaction> &&in);
  //! Add a magnetic interaction to the model (copies object)
  //! \param in magnetic interaction being added to the model builder.
  void addInteraction(const Interaction &in);
  //! Update interaction parameter
  //! \param name Name of the interaction to update.
  //! \param value New value assigned to interaction
  void updateInteraction(const std::string &name, double value);
  //! Get the classical energy associated with the model
  //! \return
  double getEnergy();
  //! Get coefficients associated with linear terms. For a ground state, each must be zero.
  //! \return array of coefficients associated with linear terms.
  Eigen::VectorXcd getFirstOrderTerms();
  //! Constructs a SpinWave object to calculate the spin-wave dispersion.
  //! \return SpinWave object used to calculate the spin-wave dispersion.
  SpinWave createElement();
private:
  Cell cell;
  InteractionsContainer interactions;
};
}
#endif
