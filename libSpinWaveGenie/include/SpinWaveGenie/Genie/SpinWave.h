#ifndef __SpinWave_H__
#define __SpinWave_H__

#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Genie/MagneticFormFactor.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Interactions/InteractionsContainer.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace SpinWaveGenie
{
#ifndef DOXYGEN_SHOULD_SKIP_THIS

struct results
{
  double weight;
  std::size_t index;
};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

//! SpinWave Class
/*!
The SpinWave Class stores the "L" matrix and calculates
the spin wave frequencies and intensities.
*/
class SPINWAVEGENIE_EXPORT SpinWave
{
public:
  //! Use SpinWaveBuilder to generate SpinWave instance
  friend class SpinWaveBuilder;
  SpinWave() = default;
  SpinWave(const Cell &cell_in, const InteractionsContainer &interactions_in);
  //! Clear dynamical matrix and any previously calculated frequencies and intensities.
  void clearMatrix();
  //! Access the cell used in this calculation.
  //! \return const reference to the cell owned by this object.
  const Cell &getCell() const;
  //! Set the location to calculate in reciprocal space.
  //! \param KX position along first axis (in rlu)
  //! \param KY position along second axis (in rlu)
  //! \param KZ position along third axis (in rlu)
  void createMatrix(double KX, double KY, double KZ);
  //! Calculate spin-wave frequencies and intensities.
  void calculate();
  //! Get calculated frequencies and intensities.
  //! \return const reference to the calculated frequencies and intensities
  const Results &getPoints() const { return VI; };

private:
  double KXP{0.}, KYP{0.}, KZP{0.};
  Cell cell;
  void calculateEigenvalues();
  void calculateWeights();
  void calculateIntensities();
  std::size_t M{0}, N{0};
  Eigen::MatrixXcd LN;
  Eigen::VectorXd SS;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  Eigen::VectorXd WW; // want to get rid of this
  Results VI;
  Eigen::MatrixXcd XY, XIN;
  InteractionsContainer interactions;
  MagneticFormFactor formFactor;
};
}
#endif /* defined(__SpinWave__) */
