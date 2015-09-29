#ifndef __SpinWave_H__
#define __SpinWave_H__

#include <iostream>
#include <utility>
#include <memory>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Genie/MagneticFormFactor.h"
#include "SpinWaveGenie/Containers/Results.h"
#include "SpinWaveGenie/Containers/Matrices.h"

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
class SpinWave
{
public:
  SpinWave() : KXP(0.0), KYP(0.0), KZP(0.0), M(0), N(0), NU(0), MI(0), IM(0), interactions{} {};
  SpinWave(const SpinWave &model);
  SpinWave &operator=(const SpinWave &other);
  //! Use SpinWaveBuilder to generate SpinWave instance
  friend class SpinWaveBuilder;
  SpinWave(Cell &cell_in, std::vector<std::unique_ptr<Interaction>> interactions_in);
  void clearMatrix();
  const Cell &getCell() const;
  void createMatrix(double KX, double KY, double KZ);
  void calculate();
  Results getPoints();

private:
  double KXP, KYP, KZP;
  Cell cell;
  void calculateEigenvalues();
  void calculateWeights();
  void calculateIntensities();
  std::size_t M, N;
  int NU, MI, IM;
  Eigen::MatrixXcd LN;
  Eigen::VectorXd SS;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
  Eigen::VectorXd WW; // want to get rid of this
  Results VI;
  Eigen::MatrixXcd XY, XIN;
  std::vector<std::unique_ptr<Interaction>> interactions;
  MagneticFormFactor formFactor;
};
}
#endif /* defined(__SpinWave__) */
