#ifndef __ExchangeInteractionSameSublattice_H__
#define __ExchangeInteractionSameSublattice_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

namespace SpinWaveGenie
{

class ExchangeInteractionSameSublattice : public Interaction
{
public:
  ExchangeInteractionSameSublattice(const std::string &name, double value, const std::string &sl_r, double min,
                                    double max);
  void updateInteraction(double value, const std::string &sl_r, double min, double max);
  void updateValue(double value_in) override;
  const std::string &getName() const override;
  void calcConstantValues(const Cell &cell) override;
  void calculateEnergy(const Cell &cell, double &energy) override;
  void calculateFirstOrderTerms(const Cell &cell, Eigen::VectorXcd &elements) override;
  void updateMatrix(const Eigen::Vector3d &K, Eigen::MatrixXcd &LN) const override;
  std::array<std::string, 2> sublattices() const override;
  std::unique_ptr<Interaction> clone() const override;

private:
  Neighbors neighbors;
  std::string name, sl_r, sl_s;
  std::size_t r, s, M;
  double value, min, max;
  std::complex<double> LNrr, LNrs;
};
}
#endif
