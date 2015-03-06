#ifndef __SpinWave_H__
#define __SpinWave_H__

#define _USE_MATH_DEFINES
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <boost/ptr_container/ptr_vector.hpp>
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
  long index;
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
  SpinWave() : KXP(0.0), KYP(0.0), KZP(0.0), M(0), N(0), NU(0), MI(0), IM(0){};
  SpinWave(const SpinWave &other)
      : KXP(other.KXP), KYP(other.KYP), KZP(other.KZP), cell(other.cell), M(other.M), N(other.N), NU(other.NU),
        MI(other.MI), IM(other.IM), LN(other.LN), SS(other.SS), ces(other.ces), WW(other.WW), VI(other.VI),
        XY(other.XY), XIN(other.XIN), interactions(other.interactions), formFactor(other.formFactor)
  {
  }
  //! Use SpinWaveBuilder to generate SpinWave instance
  friend class SpinWaveBuilder;
  SpinWave(Cell &cell_in, boost::ptr_vector<Interaction> interactions_in);
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
  boost::ptr_vector<Interaction> interactions;
  MagneticFormFactor formFactor;
};
}
#endif /* defined(__SpinWave__) */
