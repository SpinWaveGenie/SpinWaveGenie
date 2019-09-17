#ifndef __ExchangeInteractionSameSublattice_H__
#define __ExchangeInteractionSameSublattice_H__

#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT ExchangeInteractionSameSublattice : public Interaction
{
public:
  ExchangeInteractionSameSublattice(const std::string &name_in, double value_in, const std::string &sl_r_in,
                                    double min_in, double max_in);
  void updateInteraction(double value_in, const std::string &sl_r_in, double min_in, double max_in);
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
  std::size_t r, M;
  double value, min, max;
  std::complex<double> LNrr, LNrs;
};
} // namespace SpinWaveGenie
#endif
