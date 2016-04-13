#ifndef __ExchangeInteraction_H__
#define __ExchangeInteraction_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

namespace SpinWaveGenie
{

class ExchangeInteraction : public Interaction
{
public:
  ExchangeInteraction(std::string name, double value, std::string sl_r, std::string sl_s, double min, double max);
  void updateInteraction(double value, std::string sl_r, std::string sl_s, double min, double max);
  void updateValue(double value_in) override;
  const std::string &getName() const override;
  void calcConstantValues(Cell &cell) override;
  void calculateEnergy(Cell &cell, double &energy) override;
  void calculateFirstOrderTerms(Cell &cell, Eigen::VectorXcd &elements) override;
  void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN) const override;
  std::array<std::string, 2> sublattices() const override;
  std::unique_ptr<Interaction> clone() const override;

private:
  Neighbors neighbors;
  std::string name, sl_r, sl_s;
  std::size_t r, s, M;
  double value, min, max;
  std::complex<double> LNrr, LNss, LNrs, LNrsM;
};
}
#endif
